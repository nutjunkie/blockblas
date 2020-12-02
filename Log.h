#ifndef LOG_H
#define LOG_H

#include <iostream>


namespace Log
{
   inline void warn(std::string const& msg, std::ostream& os = std::cerr) 
   {
      os << "Warning: " << msg << std::endl;
   }


   inline void error(std::string const& msg, std::ostream& os = std::cerr) 
   {
      os << "ERROR: " << msg << std::endl;
   }


   inline void debug(std::string const& msg, std::ostream& os = std::cerr) 
   {
      os << "DEBUG: " << msg << std::endl;
   }

}

#endif
