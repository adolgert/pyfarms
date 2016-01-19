'''
I have scripts that do lots of little things, so I call them with
switches to tell them which little thing to do.

  python plot_obs.py --multiples
  python plot_obs.py --together

This augments the argument parser to make it a little easier to
add boolean switches that determine what the script should do.
This parser also gives some default behavior for how verbose
and quiet switches affect logging.
'''
import sys
from argparse import ArgumentParser
import logging
import logging.handlers
import unittest

class DefaultArgumentParser(ArgumentParser):
    '''This adds some default behaviors to the argument parsing.
    Default logging level is INFO.
    -v or --verbose sets logging level to DEBUG.
    -q or --quiet sets logging level to ERROR.
    --test will run tests scripts using python's unittest.
    '''

    def __init__(self,*args,**kwargs):
        '''
        Typical arguments
        \param description a text string saying what the program does. Optional.
        \param suite a function that returns a unittest test suite. Optional.
        '''
        self.tests=None
        if 'suite' in kwargs:
            self.tests=kwargs['suite']
            del kwargs['suite']
        ArgumentParser.__init__(self,*args,**kwargs)
        self.add_argument('-v','--verbose',dest='verbose',action='store_true',
                          default=False,help='print debug messages')
        self.add_argument('-q','--quiet',dest='quiet',action='store_true',
                          default=False,help='print only exceptions')
        self._functions=list()
        if self.tests:
            self.add_function('test','run unit tests in this module')



    def add_function(self,name,msg):
        '''A shorthand for created a boolean flag.'''
        self.add_argument('--%s'%name,dest=name,action='store_true',
                        default=False,help=msg)
        self._functions.append(name)


    def any_function(self):
        '''
        If none of the function switches were called, you probably
        want to show the user parser.print_help().
        '''
        did=0
        for name in self._functions:
            if getattr(self._args,name):
                did+=1
        return did>0


    def parse_args(self):
        args=ArgumentParser.parse_args(self)
        log_level=logging.INFO
        if args.verbose:
            log_level=logging.DEBUG
        elif args.quiet:
            log_level=logging.ERROR
        root_logger=logging.getLogger()
        root_logger.setLevel(log_level)

        console_handler=logging.StreamHandler(stream=sys.stdout)
        console_handler.setLevel(log_level)
        formatter=logging.Formatter(fmt="{module}.{funcName}:{lineno} {message}", datefmt="", style="{")
        console_handler.setFormatter(formatter)
        root_logger.addHandler(console_handler)

        if 'test' in dir(args) and args.test:
            unittest.TextTestRunner(verbosity=2).run(self.tests())
            sys.exit(0)
        self._args=args
        return args
