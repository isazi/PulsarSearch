
#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

// Types for the data
typedef float dataType;
const std::string dataName("float");

// DEBUG mode, prints to screen some useful information
const bool DEBUG = false;

// Memory padding
std::map< std::string, unsigned int > padding;

// Vector unit width
std::map< std::string, unsigned int > vectorWidth;

// Tuned parameters
std::map< std::string, std::map< unsigned int, std::vector< unsigned int > > > dedispersionParameters;
std::map< std::string, std::map< unsigned int, unsigned int > > transposeParameters;
std::map< std::string, std::map< unsigned int, std::map< unsigned int, std::vector< unsigned int > > > > foldingParameters;
std::map< std::string, std::map< unsigned int, std::map< unsigned int, std::vector< unsigned int > > > > snrParameters;

#endif // CONFIGURATION_HPP

