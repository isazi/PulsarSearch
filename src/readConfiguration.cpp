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
using std::getline;
#include <fstream>
using std::ifstream;
#include <map>
using std::make_pair;
#include <vector>
using std::vector;

#include <readConfiguration.hpp>
#include <utils.hpp>
using isa::utils::castToType;


void readPadding(map< string, unsigned int > & padding) {
	string temp;
	ifstream paddingFile(paddingFilename);

	while ( ! paddingFile.eof() ) {
		unsigned int middle = 0;

		getline(paddingFile, temp);
		if ( temp[0] == '#' ) {
			continue;
		}
		middle = temp.find(" ");
		padding.insert(make_pair(temp.substr(0, middle), castToType< string, unsigned int >(temp.substr(middle + 1))));
	}
}

void readVectorWidth(map< string, unsigned int > & vectorWidth) {
	string temp;
	ifstream vectorFile(vectorFilename);

	while ( ! vectorFile.eof() ) {
		unsigned int middle = 0;

		getline(vectorFile, temp);
		if ( temp[0] == '#' ) {
			continue;
		}
		middle = temp.find(" ");
		vectorWidth.insert(make_pair(temp.substr(0, middle), castToType< string, unsigned int >(temp.substr(middle + 1))));
	}
}

void readDedispersion(std::map< std::string, std::map< unsigned int, std::vector< unsigned int > > > & dedispersionParameters) {
	string temp;
	ifstream dedispersionFile(dedispersionFilename);

	while ( ! dedispersionFile.eof() ) {
		unsigned int splitPoint = 0;

		getline(dedispersionFile, temp);
		if ( temp[0] == '#' ) {
			continue;
		}
		string deviceName;
		unsigned int nrDMs = 0;
		vector< unsigned int > parameters(4);

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[0] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[1] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[2] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameters[3] = castToType< string, unsigned int >(temp);

		dedispersionParameters.insert(make_pair(deviceName, make_pair(nrDMs, parameters)));
	}
}

void readTranspose(map< string, map< unsigned int, unsigned int > > transposeParameters) {
	string temp;
	ifstream transposeFile(transposeFilename);

	while ( ! transposeFile.eof() ) {
		unsigned int splitPoint = 0;

		getline(transposeFile, temp);
		if ( temp[0] == '#' ) {
			continue;
		}
		string deviceName;
		unsigned int nrDMs = 0;
		unsigned int parameter = 0;

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameter = castToType< string, unsigned int >(temp);

		transposeParameters.insert(make_pair(deviceName, make_pair(nrDMs, parameter)));
	}
}

void readFolding(std::map< std::string, std::map< unsigned int, std::map< unsigned int, std::vector< unsigned int > > > > & foldingParameters) {
	string temp;
	ifstream foldingFile(foldingFilename);

	while ( ! foldingFile.eof() ) {
		unsigned int splitPoint = 0;

		getline(foldingFile, temp);
		if ( temp[0] == '#' ) {
			continue;
		}
		string deviceName;
		unsigned int nrDMs = 0;
		unsigned int nrPeriods = 0;
		vector< unsigned int > parameters(6);

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrPeriods = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[0] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[1] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[2] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[3] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[4] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameters[5] = castToType< string, unsigned int >(temp);

		foldingParameters.insert(make_pair(deviceName, make_pair(nrDMs, make_pair(nrPeriods, parameters))));
	}
}

void readSNR(std::map< std::string, std::map< unsigned int, std::map< unsigned int, std::vector< unsigned int > > > > & snrParameters, const string snrFilename) {
	string temp;
	ifstream snrFile(snrFilename);

	while ( ! snrFile.eof() ) {
		unsigned int splitPoint = 0;

		getline(snrFile, temp);
		if ( temp[0] == '#' ) {
			continue;
		}
		string deviceName;
		unsigned int nrDMs = 0;
		unsigned int nrPeriods = 0;
		vector< unsigned int > parameters(2);

		splitPoint = temp.find(" ");
		deviceName = temp.substr(0, splitPoint);
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrDMs = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		nrPeriods = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		splitPoint = temp.find(" ");
		parameters[0] = castToType< string, unsigned int >(temp.substr(0, splitPoint));
		temp = temp.substr(splitPoint + 1);
		parameters[1] = castToType< string, unsigned int >(temp);

		snrParameters.insert(make_pair(deviceName, make_pair(nrDMs, make_pair(nrPeriods, parameters))));
	}
}
