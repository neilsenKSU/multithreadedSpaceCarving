//
// clock.c
//
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <pthread.h>
#include "clock.h"

/* Routines for using cycle counter */

// Tell gcc not to optimize this code! 
// Otherwise, the optimizing compiler will over-optimize and mhz() will return 0
#pragma GCC optimize(0)

/* Keep track of most recent reading of cycle counter */
static uint64_t start_count = 0;

void start_counter()
{
  /* Get cycle counter */
  start_count = access_counter();
}

uint64_t access_counter()
{
  uint64_t ncycles;
  asm("rdtsc; shl $32,%%rdx; or %%rdx,%%rax"
      : "=a" (ncycles)
      : /* No input */ 
      : "%rdx");
  return(ncycles);
}

uint64_t get_counter()
{
  uint64_t end_count = access_counter();
  return (end_count - start_count);
}

double mhz()
{
  double speed;
  uint64_t cycles;

  start_counter();
  sleep(2); // sleep for 2 seconds
  cycles = get_counter();

  speed = ((double) (cycles/2) / 1.0e6);

  return(speed);
}
