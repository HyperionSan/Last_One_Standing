get_components.pl

TODO:

1)  Finish revising the authentication procedure for the various checkout tools. The script will store a list of usernames and repositories in ~/.rcl/users, and each section will determine authentication by using 
!AUTH_URL instead of !URL. If both are used, !AUTH_URL will override !URL so that would not be particularly useful right now. This will make the thornlists much more portable as users will no longer have to edit them to change the username.

2)  Confirm full Windows compatibility. Currently the script works correctly within the cygwin environment; however, the symlinks are not recognized by the Windows OS (note: cygwin does recognize them), and the logfile is seen as a single line by the Windows OS (also not an issue in cygwin). The question is do these issues represent a problem, or should the thorns only be used within the cygwin environment?

DONE 3)  Find a general solution for cvs repositories that end up in a directory like arrangments, currently the script has a special exception for anything going into 'arrangements,' but the folder could have any name and still cause issues. I think there may be a way to use the $1/$2 variables to resolve this, but I'm not sure.

DONE 4)  Add a hook to ask the user if they want to update any thorns.

DONE 5)  Add a timer functionality for checkout/update.

6)  Revise language specification, preferably into a BNL-type form.

7)  Test on all TeraGrid resources

