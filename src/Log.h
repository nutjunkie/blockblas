#ifndef LOG_H
#define LOG_H

#include <iostream>


namespace Log
{
   void warn(std::string const& msg, std::ostream& os = std::cerr);

   void error(std::string const& msg, std::ostream& os = std::cerr);

   void debug(std::string const& msg, std::ostream& os = std::cerr);
}

#endif
