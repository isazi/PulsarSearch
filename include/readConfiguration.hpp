//
// Copyright (C) 2013
// Alessio Sclocco <a.sclocco@vu.nl>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <string>
using std::string;
#include <map>
using std::map;
#include <vector>
using std::vector;


#ifndef READ_CONFIGURATION_HPP
#define READ_CONFIGURATION_HPP

void readPadding(map< string, unsigned int > & padding, const string paddingFilename);
void readVectorWidth(map< string, unsigned int > & vectorWidth, const string vectorFilename);
void readDedispersion(map< string, map< unsigned int, vector< unsigned int > > > & dedispersionParameters, const string dedispersionFilename);
void readTranspose(map< string, map< unsigned int, unsigned int > > & transposeParameters, const string transposeFilename);
void readFolding(std::map< std::string, std::map< unsigned int, std::map< unsigned int, std::vector< unsigned int > > > > & foldingParameters, const string foldingFilename);
void readSNR(std::map< std::string, std::map< unsigned int, std::map< unsigned int, std::vector< unsigned int > > > > & snrParameters, const string snrFilename);


#endif // READ_CONFIGURATION_HPP
