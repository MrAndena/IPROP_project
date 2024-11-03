#ifndef CONFIG_READER_HPP
#define CONFIG_READER_HPP

#include "data_struct.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <json.hpp>

using json = nlohmann::json;


void reader(data_struct& data);

#endif //CONFIG_READER_HPP
