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

#' Generates the the startup message
startupMessage <- function() {
	logo <- 
" _______ __    ____         __  ______ _   ______
/_  __(_) /__ / __/__ ___ _/  |/  / _ | | / / __/
 / / / / / -_)\\ \\/ -_) _ `/ /|_/ / __ | |/ / _/  
/_/ /_/_/\\__/___/\\__/\\_, /_/  /_/_/ |_|___/___/  
                      /_/               "
	return(paste0(
		logo,"v",packageVersion("tileseqMave"),"\n"
	))
}

#' Show message at startup
.onAttach <- function(lib, pkg)
{
  packageStartupMessage(startupMessage())      
  invisible()
}