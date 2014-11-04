/*!
 *  @file chrono.hpp
 *  @brief This class implements a simple chronometer.
 */

#ifndef CHRONO_HPP_
#define CHRONO_HPP_

#include <chrono>

namespace FVCode3D
{

//! Class that implements a chronometer
/*!
 * @class Chrono
 * This class implements a simple chronometer useful to measure the time.
 */
class Chrono{
public:

    //! Default constructor
    Chrono() = default;

    //! Default copy constructor
    Chrono(const Chrono &) = default;

    //! Default destructor
    ~Chrono() = default;

    //! Start the chrono
    inline void start() { M_start = Clock::now(); };

    //! Stop the chrono
    inline void stop() { M_end = Clock::now(); };

    //! Get the passed time, until this moment
    /*!
     * @return the seconds passed from the start
     * @pre call start()
     */
    inline double partial() const { return std::chrono::duration_cast<chrono_milli>(Clock::now() - M_start).count() / 1000.; }

    //! Get the passed time
    /*!
     * @return the second passed from the start to the stop
     * @pre call start()
     * @pre call stop()
     */
    inline double time() const { return std::chrono::duration_cast<chrono_milli>(M_end - M_start).count() / 1000.; }

private:
    //! Type for std::chrono::high_resolution_clock
    typedef std::chrono::high_resolution_clock Clock;

    //! Type for std::chrono::milliseconds
    typedef std::chrono::milliseconds chrono_milli;

    //! Type for std::chrono::seconds
    typedef std::chrono::seconds chrono_sec;

    //! Time at start and stop
    Clock::time_point M_start, M_end;
};

} // namespace FVCode3D

#endif /* CHRONO_HPP_ */
