
// Types for the data
typedef float dataType;
std::string dataName("float");

// DEBUG mode, prints to screen some useful information
const bool DEBUG = false;
// Memory padding
const std::map< std::string, unsigned int > padding = {
	{"GTXTitan", 32},
	{"HD7970", 64},
	{"K20", 32},
	{"Phi", 16}
};
// Vector unit width
const std::map< std::string, unsigned int > vectorWidth = {
	{"GTXTitan", 32},
	{"HD7970", 64},
	{"K20", 32},
	{"Phi", 16}
};

// Tuned parameters (APERTIF)
const std::map< std::string, const std::map< unsigned int, const std::vector< const unsigned int > > > dedispersionParameters = {
	{"GTXTitan",
		{2048,
			{32, 16, 25, 4}
		}

	},
	{"HD7970",
		{2048,
			{32, 8, 5, 4}
		}

	},
	{"K20", 
		{2048,
			{32, 16, 25, 4}
		}
	},
	{"Phi",
		{2048,
			{16, 1, 10, 4}
		}
	}
};

const std::map< std::string, const std::map< unsigned int, unsigned int > > transposeParameters = {
	{"GTXTitan",
		{2048, 32}

	},
	{"HD7970",
		{2048, 32}

	},
	{"K20", 
		{2048, 32}
	},
	{"Phi",
		{2048, 16}
	}
};

const std::map< std::string, const std::map< unsigned int, const std::map< unsigned int, std::vector< const unsigned int > > > > foldingParameters = {
	{"GTXTitan",
		{2048,
			{512, {512, 1, 1, 1, 1, 2}}
		}

	},
	{"HD7970",
		{2048,
			{512, {64, 2, 2, 1, 1, 2}}
		}

	},
	{"K20", 
		{2048,
			{512, {32, 8, 4, 1, 2, 1}}
		}
	},
	{"Phi",
		{2048,
			{512, {16, 2, 8, 4, 1, 1}}
		}
	}
};

const std::map< std::string, const std::map< unsigned int, const std::map< unsigned int, std::vector< const unsigned int > > > > snrParameters = {
	{"GTXTitan",
		{2048,
			{512, {1024, 1}}
		}

	},
	{"HD7970",
		{2048,
			{512, {128, 1}}
		}

	},
	{"K20", 
		{2048,
			{512, {32, 16}}
		}
	},
	{"Phi",
		{2048,
			{512, {16, 2}}
		}
	}
};