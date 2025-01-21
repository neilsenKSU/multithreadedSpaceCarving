//
// clock.h - Routines for using cycle counter
//

/* Start the counter, set start_count */
void start_counter();

/* Access cycle counter */
uint64_t access_counter();

/* Get # elapsed cycles since counter started */
uint64_t get_counter();

/* Compute machine speed in MHz; e.g., 2600 MHz = 2.6 GHz */
double mhz();
