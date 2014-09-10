#!/usr/bin/env python

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# Create a Python logger with the right parameters.

import sys
import logging

logLevel = logging.DEBUG
logFormat = '[%(levelname)s] %(asctime)-15s [%(module)s %(lineno)d] %(message)s'
logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
logger = logging.getLogger ('MatchAnnot')

# Copyright (C) 2011 Genome Research Limited
#
# This library is free software. You can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
