#include <iostream>
#include "Log.h"


namespace Log
{
   void warn(std::string const& msg, std::ostream& os) 
   {
      os << "Warning: " << msg << std::endl;
   }


   void error(std::string const& msg, std::ostream& os) 
   {
      os << "ERROR: " << msg << std::endl;
   }


   void debug(std::string const& msg, std::ostream& os) 
   {
      os << "DEBUG: " << msg << std::endl;
   }
}
