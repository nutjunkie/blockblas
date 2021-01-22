#ifndef TIMER_H
#define TIMER_H
#include <chrono>
#include <string>
#include <iomanip>
#include <sstream>


class Timer {

   public:
      typedef std::chrono::steady_clock clock;

      Timer() : m_elapsed(0.0) { }

      void start() 
      {
         m_elapsed = 0.0;
         m_start = clock::now();
      }

      void restart() 
      {
         m_start = clock::now();
      }

      double stop() 
      {
         std::chrono::time_point<clock> stop = clock::now();
         std::chrono::milliseconds ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(stop-m_start);
         m_elapsed += 0.001*ms.count();
         return m_elapsed;
      }

      std::string format(double const time) const
      {
         std::stringstream s;         
         s << std::fixed << std::showpoint << std::setprecision(3);
         s << time << " s";
         return s.str();
      }

      std::string format() const 
      {
         return format(m_elapsed);
      }

   private:
      double m_elapsed;
      std::chrono::time_point<clock> m_start;
};

#endif
