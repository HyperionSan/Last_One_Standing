#!/usr/bin/env python

# this script will repeatedly call GetComponents and SimFactory using different timestamps for
# the checkout, in attempt to isolate when an error was introduced into the source tree.

import os
import sys
from datetime import date, timedelta
from optparse import OptionParser
from shutil import rmtree

checkout_dates = []

orig_dir = os.getcwd()

def main():
    usage = "usage: %prog [options] thornlist"
    parser = OptionParser(usage=usage)
    parser.add_option("-r", "--root", action="store",
                        dest="root", default="Cactus",
                        help="set root for Cactus checkout")
    parser.add_option("-t", "--test-script", action="store",
                        dest="test", help="provide custom script to build and test Cactus")
    parser.add_option("-n", "--number", action="store",
                        dest="number", default=5, help="how many copies to keep. defaults to 5")

    (options, args) = parser.parse_args()
    thornlist = os.path.join(orig_dir, args[0])

    d = date.today()
    t = timedelta(days=30)
    done = False
    passed = False

    while not done:
        checkout_source(options, thornlist, d)
        passed = test_source(options, thornlist, d)
        (d, t, done) = change_date(d, t, passed)

    print "Successful compilation at %s" % d.isoformat()


def checkout_source(options, thornlist, d):
    cmd = './GetComponents -a --root=%s.%s --date %s %s' % (options.root, d.isoformat(), d.isoformat(), thornlist)
    err = os.system(cmd)
    if err != 0:
        print "Error: Could not complete checkout from %s" % d.isoformat()
        sys.exit(1)
    checkout_dates.append(d)
    # remove extra checkouts
    if len(checkout_dates) > options.number:
        rmtree(options.root+"."+checkout_dates[0].isoformat())

def test_source(options, thornlist, d):
    # default if no test script is supplied
    if not options.test:
        os.chdir(options.root+'.'+d.isoformat())
        os.system('cp simfactory/udb.example.pm simfactory/udb.pm')
        err = os.system('./simfactory/sim build --thornlist=%s' % thornlist)
        if err == 0:
            passed = True
        else:
            passed = False
        os.chdir(orig_dir)

    # if test script supplied, run it
    else:
        err = os.system(options.test)
        if err == 0:
            passed = True
        else:
            passed = False

    return passed

def change_date(d, t, passed):
    if passed is False:
        # only decrement the timedelta if we have already found a good checkout date
        # at this point the timedelta will not have its starting value
        if t != timedelta(days=30):
            t /= 2
        d -= t
    else:
        # first make sure we are not at today
        if d == date.today():
            done = True
            #print "I'm DONE!!"
        # otherwise decrement the timedelta and modify Date
        else:
            #print "I'm NOT DONE!!"

            t /= 2
            d += t
            # compare with latest checkout date
            if d == checkout_dates[-1]:
                done = True

    return d, t, done

######################################################
if __name__ == '__main__':
    main()