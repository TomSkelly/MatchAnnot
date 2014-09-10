#!/usr/bin/env python

# Keep track of the best value we're seen so far in a series of calls.

# This is pretty elaborate. But this way we're got all the gubbins in
# one place, and the calls to it in the main body of the program are
# simple. The parent script was starting to fill up with "find the
# best one" clutter.

# One behavioral nuance: If the new candidate's value is *equal* to
# the current best, we *do not* accept it as a new best.

import os
import sys

class Best (object):
    '''Keep track of the best <whatever> we're seen.'''

    def __init__ (self, reverse=False):
        self.value   = None
        self.which   = None
        self.reverse = reverse

    def update (self, value, which):           # returns True if we have a new winner

        if self.value is None:
            self.value = value
            self.which = which
            return True

        if not self.reverse:                   # if big is good
            if value > self.value:
                self.value = value
                self.which = which
                return True

        else:
            if value < self.value:
                self.value = value
                self.which = which
                return True

        return False

