 /*@@
   @file      Sockets.c
   @date      Wed Sep 13 20:39:15 2000
   @author    Tom Goodale
   @desc
              Routines which deal with sockets.
              These should probably move into thorn Socket at some point if
              they aren't already there.
   @enddesc
   @version   $Id: Sockets.c 101 2011-08-17 14:55:38Z knarf $
 @@*/

#include "cctk.h"
#include "util_Network.h"
#include "CactusBase/IOUtil/src/ioutil_AdvertisedFiles.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif /* HAVE_SYS_TIME_H */
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif /* HAVE_SYS_TYPES_H */
#ifdef HAVE_SYS_SOCKET_H
#include <sys/socket.h>
#endif /* HAVE_SYS_SOCKET_H */
#ifdef HAVE_NETINET_IN_H
#include <netinet/in.h>
#endif /* HAVE_NETINET_IN_H */
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif /* HAVE_NETDB_H */
#ifdef HAVE_ARPA_INET_H
#include <arpa/inet.h>
#endif /* ARPA_INET_H */
#ifdef HAVE_WINSOCK2_H
#include <winsock2.h>
#endif /* HAVE_WINSOCK2_H */
#include <errno.h>

#ifndef WIN32
#include <signal.h>
#endif

#include "IsoSurfacerInit.h"
#include "IsoSurfacerGH.h"
#include "util_String.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusPUGHIO_IsoSurfacer_Sockets_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

#ifndef HAVE_SOCKET_TYPE
#define SOCKET int
#endif

#ifdef SOCKET_ERROR
#define ERROR_CHECK(a)  ((a) == SOCKET_ERROR)
#else
#define ERROR_CHECK(a)  ((a) < 0)
#endif

#ifdef HAVE_WINSOCK2_H
#define CLOSESOCKET(a) closesocket(a)
#else
#define CLOSESOCKET(a) close(a)
#endif

#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0
#endif


typedef enum {closed, open} isoSocketState;

typedef struct ISOSocket
{
  struct ISOSocket *prev;
  struct ISOSocket *next;
  unsigned long filedes;
  isoSocketState state;
} isoSocket;

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int byteswap(void *buf,CCTK_INT8 nelements,int elementsize);

IsoCommand *Iso_PollCommand(const cGH *GH,IsoCommand *cmd);
int Iso_Write(isoSocket *connection, const char *buffer, size_t count);

static SOCKET Iso_MakeSocket (unsigned long port, int *hunt);

static isoSocket *SocketCreate(unsigned long int filedes, isoSocket **list);
static void SocketDestroy(isoSocket *this, isoSocket **list);
static void SocketClose(isoSocket *this);
int Iso_SetupServer(const cGH *GH,isosurfacerGH *myGH,
                    int dataport, int controlport, int queue_size, int hunt);
void Iso_ShutdownServer(CCTK_ARGUMENTS);
int Iso_Poll(const cGH *GH, long sec, long usec);
int Iso_Read(isoSocket *connection, char *buffer, size_t count);
int IsoWriteDataToClients(const char *metadata,
                          CCTK_INT8 size,
                          IsoType type,
                          void *dataP);

static int InitialiseTCP(void);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/
/* OK, maybe this stuff should go into the isosurfacer GH?
   Just general fear of global variables... putting it into the
   GH uses the same amount of space as the static vars, but gives
   us a sort of a "namespace".
*/
/* Active file descriptors */
static fd_set active_fd_set;

/* Main server socket */
static SOCKET datasock;
static SOCKET controlsock;

static SOCKET minsock;
static SOCKET maxsock;

static isoSocket *controlsocklist = NULL;
static isoSocket *datasocklist = NULL;

static int chosen_controlport = 0;
static int chosen_dataport = 0;

#define FILENAME_TEMPLATE "fileXXXXXX"
static char *advertised_filename;


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/


 /*@@
   @routine    Iso_SetupServer
   @date       Wed Sep 13 20:39:15 2000
   @author     Tom Goodale
   @desc
   Creates a socket to listen on.
   @enddesc
@@*/
int Iso_SetupServer(const cGH *GH,isosurfacerGH *myGH,int dataport, int controlport, int queue_size, int hunt)
{
  int realdport;
  int realcport;
  ioAdvertisedFileDesc advertised_file;
  FILE *advertised_file_fd;
  char hostname[1025];

  advertised_filename = Util_Strdup(FILENAME_TEMPLATE);
  if (advertised_filename == NULL)
    CCTK_VWarn (0, __LINE__, __FILE__, "IsoSurfacer", "Cannot allocate memory foradvertised file name");
  myGH = myGH;

  if(CCTK_MyProc(GH) != 0) return 0; /* not the root processor */

  /* Some systems need special logic for starting up TCP. */
  InitialiseTCP();

  /* Create the data socket and set it up to accept connections. */
  datasock = Iso_MakeSocket (dataport, hunt ? &realdport : NULL);

  if (ERROR_CHECK(listen (datasock, queue_size)))
  {
    perror ("listen");
    exit (EXIT_FAILURE);
  }

  /* Create the control socket and set it up to accept connections. */
  controlsock = Iso_MakeSocket (controlport, hunt ? &realcport : NULL);

  if (ERROR_CHECK(listen (controlsock, queue_size)))
  {
    perror ("listen");
    exit (EXIT_FAILURE);
  }

  chosen_controlport = hunt ? realcport:controlport;
  chosen_dataport = hunt ? realdport:dataport;
  minsock = datasock > controlsock ? controlsock : datasock;
  maxsock = datasock > controlsock ? datasock    : controlsock;

  /* Initialize the set of active sockets. */
  FD_ZERO (&active_fd_set);
  FD_SET (controlsock, &active_fd_set);
  FD_SET (datasock, &active_fd_set);

  Util_GetHostName(hostname, sizeof (hostname));
  CCTK_VInfo (CCTK_THORNSTRING,
              "Isosurfacer listening for connections\n"
              "                   host '%s' control port %d data port %d",
              hostname, chosen_controlport, chosen_dataport);

#ifdef HAVE_MKSTEMP
  advertised_file_fd = fdopen (mkstemp (advertised_filename), "w");
#else
  advertised_file_fd = tmpnam (advertised_filename) ?
                       fopen (advertised_filename, "w") : NULL;
#endif
  if (advertised_file_fd)
  {
    fprintf (advertised_file_fd, "Hostname:     %s\n", hostname);
    fprintf (advertised_file_fd, "Control port: %d\n", chosen_controlport);
    fprintf (advertised_file_fd, "Data port:    %d\n", chosen_dataport);
    fclose (advertised_file_fd);

    /* FIXME: this can go after the old filename scheme has gone */
    advertised_file.slice = "";
    advertised_file.thorn = CCTK_THORNSTRING;
    advertised_file.varname = "All variables";
    advertised_file.description = "Streamed isosurface geometry data";
    advertised_file.mimetype = "data/x-streamed-isosurfaces";

    IOUtil_AdvertiseFile (GH, advertised_filename, &advertised_file);
  }
  else
  {
    CCTK_WARN (2, "Couldn't create unique temporary file ! "
                  "Isosurfacer data streaming was not advertised.");
  }

  return (0);
}

 /*@@
   @routine    Iso_ShutdownServer
   @date       Wed Sep 13 20:39:15 2000
   @author     Tom Goodale
   @desc
   Closes all sockets we are interested in.
   @enddesc
@@*/
void Iso_ShutdownServer(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  SOCKET i;

  /* Close all sockets in our active set */
  for(i = maxsock; i >= minsock; i--)
  {
    if(FD_ISSET(i, &active_fd_set))
    {
      CLOSESOCKET(i);
    }
  }

  /* remove advertised file and free filename */
  if (strcmp (advertised_filename, FILENAME_TEMPLATE))
  {
    remove (advertised_filename);
  }
  free (advertised_filename);
  return;
}

 /*@@
   @routine    Iso_Poll
   @date       Wed Sep 13 20:39:15 2000
   @author     Tom Goodale
   @desc
   Main workhorse routine.
   Looks for activity on any of the sockets which are open
   and dispatches work.
   @enddesc
@@*/
int Iso_Poll(const cGH *GH, long sec, long usec)
{

#ifdef HAVE_SOCKLEN_T
  socklen_t size;
#else
  int size;
#endif

  fd_set read_fd_set;

  struct sockaddr_in clientname;
  struct timeval timeout;
  struct timeval *real_timeout;


  if(CCTK_MyProc(GH) != 0) return 0; /* not the root processor */
  if(sec >=0)
  {
    timeout.tv_sec = sec;
    timeout.tv_usec = usec;

    real_timeout = &timeout;
  }
  else
  {
    real_timeout = NULL;
  }

  /* Check if any input is available on one or more active sockets. */
  read_fd_set = active_fd_set;

  if (ERROR_CHECK(select (FD_SETSIZE, &read_fd_set, NULL, NULL, real_timeout)))
  {
    perror ("select");
    abort ();
  }

  /* Service all the sockets with input pending. */
  if (FD_ISSET (datasock, &read_fd_set))
  {
    /* Connection request on data socket. */
    SOCKET new;
    size = sizeof (clientname);
    new = accept (datasock,
                  (struct sockaddr *) &clientname,
                  &size);
    if (ERROR_CHECK(new))
    {
      perror ("accept");
      abort ();
    }
    fprintf (stderr,
             "Isosurfacer: data connect to port %u from host %s, port %hd.\n",
             chosen_dataport,
             inet_ntoa (clientname.sin_addr),
             ntohs (clientname.sin_port));

    SocketCreate(new, &datasocklist);
  }

  if (FD_ISSET (controlsock, &read_fd_set))
  {
    /* Connection request on data socket. */
    SOCKET new;
    size = sizeof (clientname);
    new = accept (controlsock,
                  (struct sockaddr *) &clientname,
                  &size);
    if (ERROR_CHECK(new))
    {
      perror ("accept");
      abort ();
    }
    fprintf (stderr,
             "Isosurfacer: control connect to port %u from host %s, port %hd.\n",
             chosen_controlport,
             inet_ntoa (clientname.sin_addr),
             ntohs (clientname.sin_port));

    SocketCreate(new, &controlsocklist);
  }

  return 0;
}

IsoCommand *Iso_PollCommand(const cGH *GH,IsoCommand *cmd)
{
  fd_set read_fd_set;
  struct timeval timeout;
  struct timeval *real_timeout;

  isoSocket *this;


  if(CCTK_MyProc(GH) != 0) return 0; /* not the root processor */

  /* always a timeout of 0 */
  timeout.tv_sec = 0;
  timeout.tv_usec = 0;
  real_timeout = &timeout;

  /* Check if any input is available on one or more active sockets. */
  read_fd_set = active_fd_set;

  if (ERROR_CHECK(select (FD_SETSIZE, &read_fd_set, NULL, NULL, real_timeout))){
    perror ("select");
    abort ();
  }
  /* get next command on the control socket list */
  for(this = controlsocklist; this; this = this->next){
    if(FD_ISSET (this->filedes, &read_fd_set)){
      int r;
      /* puts("StartIsoRead"); */
      r=Iso_Read(this,cmd->buffer,64+64+64);
      /* puts("EndIsoRead"); */
      /* Data arriving on an already-connected socket. */
      if (r <= 0){
#if 0
        puts("Destroy Socket");
#endif
        SocketDestroy(this, &controlsocklist);
#if 0
        puts("done destroy");
        if(controlsocklist==NULL) puts("*********It did the right thing");
#endif
        if(controlsocklist) puts("*******Failed to destroy socket*******");
        return NULL;
      }
      else return cmd; /* return immediately with the new command (first available) */
    }
  }
  return NULL;
}

int IsoWriteDataToClients(const char *metadata,
                          CCTK_INT8 size,
                          IsoType type,
                          void *dataP)
{
  CCTK_INT4 *data=(CCTK_INT4*)dataP;
  int retval;
  CCTK_INT4 datatype=type;
  CCTK_INT8 datasize=size;
  isoSocket *this;
  isoSocket *next;


  /* Convert the data to network byte order */
  byteswap(data,size,4);
  /* for(i = 0; i < size; i++)
   {
     I wish htons/htonl actually worked as you'd expect them to
     data[i] = htons(data[i]);
   }*/
  /*  printf("Pre-Swap Datasize=%lld:%llx\n",datasize,datasize); */
  byteswap(&datatype,1,4);
  byteswap(&datasize,1,8);
  /* printf("PostSwap Datasize=%lld:%llx\n",datasize,datasize); */
  /*  unfortunately short is 16 bits and long is 32 bits on some systems
  datatype = htons(datatype);
  datasize = htonl(datasize); */

  /* Send data down all connected sockets. */
  for(this = datasocklist; this; this = next)
  {
    next = this->next;

    retval=Iso_Write(this, (void*)(&datasize), 8);
    /* printf("RetVal on size write=%d\n",retval); */
    Iso_Write(this, (void*)(&datatype), 4);
    Iso_Write(this, metadata, 64);
    retval = Iso_Write(this, (void *)data, size*4);

    if(retval < 0)
    {
      SocketDestroy(this, &datasocklist);
    }
  }

  /* Convert the data back */
  byteswap(data,size,4);
#if 0
  for(i = 0; i < size; i++)
  {
    data[i] = htons(data[i]);
  }
#endif

  return (0);
}


 /*@@
   @routine    Iso_Write
   @date       Fri Sep 15 18:47:41 2000
   @author     Tom Goodale
   @desc
   Writes part or all of an Iso data connection.
   @enddesc
@@*/
int Iso_Write(isoSocket *connection, const char *buffer, size_t count)
{
  int retval;
  size_t bytes_sent;
  size_t halfcount;
  int done;
  int tmp;
  fd_set this_set;

  if(connection->state == open)
  {
    FD_ZERO (&this_set);
    done = 0;
    bytes_sent = 0;
    while(!done)
    {
      /* Try and send the data.  Make sure we don't get a SIGPIPE if the
       * other end is closed.
       */
      retval = send(connection->filedes,
                    buffer+bytes_sent,
                    count-bytes_sent,
                    MSG_NOSIGNAL);

      /* Check for errors */
      if(ERROR_CHECK(retval))
      {
        switch(errno)
        {
#ifdef EMSGSIZE
          case EMSGSIZE : /* The socket requires that message be sent atomically,
                           *  and the size of the message to be sent made this
                           * impossible.
                           */
            /* Split message in two and try again. */
            halfcount = count/2;
            retval = Iso_Write(connection, buffer, halfcount);
            if(retval > -1)
            {
              tmp = Iso_Write(connection, buffer+halfcount, count-halfcount);
              if(tmp > -1)
              {
                retval += tmp;
              }
              else
              {
                retval = tmp;
              }
            }
            break;
#endif

#ifdef ENOBUFS
          case ENOBUFS :
            /* The output queue for a network interface was full.  This generally indicates that the  interface  has
             * stopped sending, but may be caused by transient congestion.  (This cannot occur in Linux, packets are
             * just silently dropped when a device queue overflows.)
             */
            retval = 0;
            break;
#endif
          case EBADF    : /* An invalid descriptor was specified. */
#ifdef ENOTSOCK
          case ENOTSOCK : /* The filedescriptor is not a socket. */
#endif
#ifdef ECONNRESET
          case ECONNRESET : /* Connection reset by peer */
#endif
          case EPIPE : /* The local end has been shut down on a connection oriented socket.  In this case the process will also
                        *  receive a SIGPIPE unless MSG_NOSIGNAL is set.
                        */
            SocketClose(connection);
            break;

          case EINTR : /* A signal occurred.*/

          case ENOMEM :  /* No memory available. */

          case EINVAL : /* Invalid argument passed. */

          case EFAULT   : /* An invalid user space address was specified for a parameter. */
            break;

#if 0
          case EAGAIN      :
          case EWOULDBLOCK :
            /* The socket is marked non-blocking and the requested operation would block. */
#endif
          default :
            fprintf(stderr, "Unhandled error '%s':  this should never happen - %s %d", strerror (errno), __FILE__, __LINE__);
        }
      }
      else
      {
        bytes_sent += retval;
      }

      if(!(ERROR_CHECK(retval)) && bytes_sent < count)
      {
        /* Wait until can write */
        FD_SET (connection->filedes, &this_set);
        if (ERROR_CHECK(select (FD_SETSIZE, NULL, &this_set, NULL, NULL)))
        {
          perror ("select");
          CCTK_Abort(NULL, EXIT_FAILURE);
        }
      }
      else
      {
        done = 1;
      }
    };
  }
  else
  {
    retval = -1;
  }

#ifdef Iso_DEBUG
  if(retval != -1)
  {
    fprintf (stderr, "Wrote: `%s'\n", buffer);
  }
  else
  {
    fprintf(stderr, "Couldn't write to the socket 8-(\n");
  }
#endif

  return retval;
}

 /*@@
   @routine    Iso_Read
   @date       Mon Sep 18 10:14:03 2000
   @author     Tom Goodale
   @desc
   Reads part or all of an Iso request.
   @enddesc
@@*/
int Iso_Read(isoSocket *connection, char *buffer, size_t count)
{
  int retval;

  if(connection->state == open)
  {
    /* Currently don't do anything fancy. */
    retval = recv(connection->filedes, buffer, count,MSG_NOSIGNAL);

#ifdef Iso_DEBUG
    fprintf (stderr, "Read: `%s'\n", buffer);
#endif
  }
  else
  {
    retval = -1;
  }

  return retval;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    byteswap
   @date
   @author     John Shalf
   @desc

   @enddesc
@@*/
static int byteswap(void *buf,CCTK_INT8 nelements,int elementsize)
{
  char *buffer;
  CCTK_INT8 i;
  int s,d;
#ifdef WORDS_BIGENDIAN
  return 0;
#else
  buffer=(char *)buf;
  if(elementsize<=1) return 0;

  for(i=0;i<nelements;i++,buffer+=elementsize){
    /* do the swap thing on each element */
    for(s=0,d=elementsize-1;s<d;s++,d--){
      char        c=buffer[s];
      buffer[s]=buffer[d];
      buffer[d]=c;
    }
  }
  return 1;
#endif
}

 /*@@
   @routine    Iso_MakeSocket
   @date       Wed Sep 13 20:39:15 2000
   @author     Tom Goodale
   @desc
   Creates a socket.
   @enddesc
@@*/
int Iso_MakeSocket (unsigned long port, int *hunt)
{
  int this_sock;
  struct sockaddr_in name;
  int opt;
  int done;
  int realport;
  int error;

  /* Create the socket. */
  this_sock = socket (PF_INET, SOCK_STREAM, 0);
  if (ERROR_CHECK(this_sock))
  {
    perror ("socket");
    CCTK_Abort(NULL, EXIT_FAILURE);
  }

  /* Try to reuse the port if possible */
#ifdef SO_REUSEADDR
  opt = 1;
  setsockopt(this_sock, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));
#endif

  done = 0;
  realport = port;

  while(!done)
  {
    /* Give the socket a name. */
    name.sin_family = AF_INET;
    name.sin_port = htons (realport);
    name.sin_addr.s_addr = htonl (INADDR_ANY);
    if (ERROR_CHECK(error = bind (this_sock, (struct sockaddr *) &name, sizeof (name))))
    {
      if(hunt)
      {
        /* Hunt for a new port */
        fprintf(stderr, "Port %u taken, trying %u\n", realport, realport+1);
        realport++;
      }
      else
      {
        done = 1;
      }
    }
    else
    {
      done = 1;
    }
  }

  if(ERROR_CHECK(error))
  {
    perror ("bind");
    CCTK_Abort(NULL,EXIT_FAILURE);
  }
  else if(hunt)
  {
    *hunt = realport;
  }

  return this_sock;
}

 /*@@
   @routine    SocketCreate
   @date       Thu Sep 21 15:03:32 2000
   @author     Tom Goodale
   @desc
   Creates an isoSocket structure and links it to the list.
   @enddesc
@@*/
static isoSocket *SocketCreate(unsigned long int filedes, isoSocket **list)
{
  isoSocket *this;

  this = malloc(sizeof (isoSocket));

  if(this)
  {
    this->filedes = filedes;
    this->state = open;
    this->prev = NULL;
    this->next = *list;

    if(*list)
    {
      (*list)->prev = this;
    }

    *list = this;

    FD_SET (this->filedes, &active_fd_set);
  }

  return this;
}

 /*@@
   @routine    SocketDestroy
   @date       Thu Sep 21 15:04:03 2000
   @author     Tom Goodale
   @desc
   Destroys an isoSocket structure and removes it from the list.
   @enddesc
@@*/
static void SocketDestroy(isoSocket *this, isoSocket **list)
{
  if(this)
  {
    if(this->state == open)
    {
      SocketClose(this);
    }

    if(this->next)
    {
      this->next->prev = this->prev;
    }

    if(this->prev)
    {
      this->prev->next=this->next;
    }
    else
    {
      *list = this->next;
    }

    free(this);
  }
}

 /*@@
   @routine    SocketClose
   @date       Thu Sep 21 15:04:27 2000
   @author     Tom Goodale
   @desc
   Closes the socket associated with an isoSocket structure.
   @enddesc
@@*/
static void SocketClose(isoSocket *this)
{
  if(this)
  {
    if(this->state == open)
    {
      CLOSESOCKET(this->filedes);
      this->state=closed;
      FD_CLR (this->filedes, &active_fd_set);
    }
  }
}

/******************************************************************************
 ******************************************************************************
 ******************************************************************************/

/* Special code for starting up the socket layer. */

#ifdef HAVE_WINSOCK2_H
#ifndef _M_IX86
#define _M_IX86 400
#endif

#include <sys/types.h>
#include <windows.h>
#include <stdio.h>

static int InitialiseTCP(void)
{
  WORD wVersionRequested;
  WSADATA wsaData;
  int errnumber;

  wVersionRequested = MAKEWORD( 2, 0 );

  errnumber  = WSAStartup( wVersionRequested, &wsaData );

  if (errnumber)
  {
    fprintf(stderr, "Couldn't start Windows socket layer\n");
  }
  else
  {
    printf("Windows socket layer initialized.\n");
  }

  return errnumber;
}

#else

static int InitialiseTCP(void)
{
  /* Make sure we ignore SIGPIPE */
  signal(SIGPIPE,SIG_IGN);
  return 0;
}

#endif /* defined HAVE_WINSOCK2_H */
