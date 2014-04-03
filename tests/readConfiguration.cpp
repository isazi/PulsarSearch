// Copyright 2014 Alessio Sclocco <a.sclocco@vu.nl>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <exception>
using std::exception;
#include <map>
#include <vector>

#include <ArgumentList.hpp>
using isa::utils::ArgumentList;
#include <configuration.hpp>
#include <readConfiguration.hpp>


int main(int argc, char *argv[]) {

  try {
    ArgumentList args(argc, argv);

    readPadding(padding, args.getSwitchArgument< string >("-padding_file"));
    readVectorWidth(vectorWidth, args.getSwitchArgument< string >("-vector_file"));
    readDedispersion(dedispersionParameters, args.getSwitchArgument< string >("-dedispersion_file"));
    readTranspose(transposeParameters, args.getSwitchArgument< string >("-transpose_file"));
    readFolding(foldingParameters, args.getSwitchArgument< string >("-folding_file"));
    readSNR(snrParameters, args.getSwitchArgument< string >("-snr_file"));
  } catch ( exception &err ) {
    cerr << err.what() << endl;
    return 1;
  }

  // Checking the padding
  cout << "Padding" << endl;
  for ( std::map< std::string, unsigned int >::const_iterator iter0 = padding.begin(); iter0 != padding.end(); ++iter0 ) {
    cout << "\t" << (*iter0).first << " " << (*iter0).second << endl;
  }
  cout << endl;
  // Checking the vector
  cout << "Vector" << endl;
  for ( std::map< std::string, unsigned int >::const_iterator iter0 = vectorWidth.begin(); iter0 != vectorWidth.end(); ++iter0 ) {
    cout << "\t" << (*iter0).first << " " << (*iter0).second << endl;
  }
  cout << endl;
  // Checking the dedispersion
  cout << "Dedispersion"<< endl;
  for ( std::map< std::string, std::map< unsigned int, std::vector< unsigned int > > >::const_iterator iter0 = dedispersionParameters.begin(); iter0 != dedispersionParameters.end(); ++iter0 ) {
    cout << "\t" << (*iter0).first << " ";

    for ( std::map< unsigned int, std::vector< unsigned int > >::const_iterator iter1 = (*iter0).second.begin(); iter1 != (*iter0).second.end(); ++iter1 ) {
      cout << (*iter1).first << " ";

      for ( std::vector< unsigned int >::const_iterator iter2 = (*iter1).second.begin(); iter2 != (*iter1).second.end(); ++iter2 ) {
        cout << *iter2 << " ";
      }
      cout << "; ";
    }
    cout << endl;
  }
  cout << endl;
  // Checking the transpose
  cout << "Transpose" << endl;
  for ( std::map< std::string, std::map< unsigned int, unsigned int > >::const_iterator iter0 = transposeParameters.begin(); iter0 != transposeParameters.end(); ++iter0 ) {
    cout << "\t" << (*iter0).first << " ";

    for ( std::map< unsigned int, unsigned int >::const_iterator iter1 = (*iter0).second.begin(); iter1 != (*iter0).second.end(); ++iter1 ) {
      cout << (*iter1).first << " " << (*iter1).second << " ;";
    }
    cout << endl;
  }
  cout << endl;
  // Checking the folding
  cout << "Folding"<< endl;
  for ( std::map< std::string, std::map< unsigned int, std::map< unsigned int, std::vector< unsigned int > > > >::const_iterator iter0 = foldingParameters.begin(); iter0 != foldingParameters.end(); ++iter0 ) {
    cout << "\t" << (*iter0).first << " ";

    for ( std::map< unsigned int, std::map< unsigned int, std::vector< unsigned int > > >::const_iterator iter1 = (*iter0).second.begin(); iter1 != (*iter0).second.end(); ++iter1 ) {
      cout << (*iter1).first << " ";

      for ( std::map< unsigned int, std::vector< unsigned int > >::const_iterator iter2 = (*iter1).second.begin(); iter2 != (*iter1).second.end(); ++iter2 ) {
        cout << (*iter2).first << " ";

        for ( std::vector< unsigned int >::const_iterator iter3 = (*iter2).second.begin(); iter3 != (*iter2).second.end(); ++iter3 ) {
          cout << *iter3 << " ";
        }
        cout << "; ";
      }
    }
    cout << endl;
  }
  cout << endl;
  // Checking the snr
  cout << "SNR"<< endl;
  for ( std::map< std::string, std::map< unsigned int, std::map< unsigned int, std::vector< unsigned int > > > >::const_iterator iter0 = snrParameters.begin(); iter0 != snrParameters.end(); ++iter0 ) {
    cout << "\t" << (*iter0).first << " ";

    for ( std::map< unsigned int, std::map< unsigned int, std::vector< unsigned int > > >::const_iterator iter1 = (*iter0).second.begin(); iter1 != (*iter0).second.end(); ++iter1 ) {
      cout << (*iter1).first << " ";

      for ( std::map< unsigned int, std::vector< unsigned int > >::const_iterator iter2 = (*iter1).second.begin(); iter2 != (*iter1).second.end(); ++iter2 ) {
        cout << (*iter2).first << " ";

        for ( std::vector< unsigned int >::const_iterator iter3 = (*iter2).second.begin(); iter3 != (*iter2).second.end(); ++iter3 ) {
          cout << *iter3 << " ";
        }
        cout << "; ";
      }
    }
    cout << endl;
  }
  cout << endl;

  return 0;
}
