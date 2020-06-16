# Copyright (C) 2018  Jochen Weile, Roth Lab
#
# This file is part of tileseqMave.
#
# tileseqMave is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tileseqMave is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with tileseqMave.  If not, see <https://www.gnu.org/licenses/>.

#' Globally registers a logger object
#' 
#' @param logger a yogilogger object
#' @return nothing
#' @export
registerLogger <- function(logger) {
	stopifnot(inherits(logger,"yogilogger"))
	options(tileseqMave.logger=logger)
}

#' convenience logging function at INFO level
#' 
#' Checks whether a logger has been registered globally and if so, uses it.
#' Otherwise writes the message to stdout.
#' 
#' @param ... any objects or messages to be logged
#' @return nothing
#' @export
logInfo <- function(...) {
	logger <- getOption("tileseqMave.logger")
	if (!is.null(logger)) {
		logger$info(...)
	} else {
		do.call(cat,c(list(...),"\n"))
	}
}

#' convenience logging function at WARN level
#' 
#' Checks whether a logger has been registered globally and if so, uses it.
#' Otherwise writes the message to stdout.
#' 
#' @param ... any objects or messages to be logged
#' @return nothing
#' @export
logWarn <- function(...) {
	logger <- getOption("tileseqMave.logger")
	if (!is.null(logger)) {
		logger$warning(...)
	} else {
		do.call(cat,c("Warning:",list(...),"\n"))
	}
}

#' Convenience function for recording the package version into the logfile.
#' @export
logVersion <- function() {
	logInfo("Running tileseqMave version",packageVersion("tileseqMave"))
}

#There is no convenience method for error level, as errors should always be thrown
# and propagated through the call stack.