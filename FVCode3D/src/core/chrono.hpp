/*!
 *	@file chrono.hpp
 *	@brief This class implements a simple chronometer.
 */

#ifndef CHRONO_HPP_
#define CHRONO_HPP_

#include <ctime>

//! Class that implements a chronometer
/*!
 * @class Chrono
 * This class implements a simple chronometer useful to measure the time.
 */
class Chrono{
public:

	//! Start the chrono
	inline void start() { M_start = clock(); };

	//! Stop the chrono
	inline void stop() { M_end = clock(); };

	//! Get the passed time, until this moment
	/*!
	 * @return the seconds passed from the start
	 * @pre call start()
	 */
	inline double partial() const { return (double) ( (clock() - M_start) / CLOCKS_PER_SEC ); };

	//! Get the passed time
	/*!
	 * @return the second passed from the start to the stop
	 * @pre call start()
	 * @pre call stop()
	 */
	inline double time() const { return (double) ( (M_end - M_start) / CLOCKS_PER_SEC ); };

private:

	//! Time at start and stop
	clock_t M_start, M_end;

};

#endif /* CHRONO_HPP_ */
