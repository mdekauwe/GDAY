/*
Computes the simple moving average of a series of numbers.
Taken from:
    http://rosettacode.org/wiki/Averages/Simple_moving_average
*/



#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "simple_moving_average.h"

sma_result sma(enum Action action, ...)
{
     va_list vl;
     sma_result r;
     sma_obj *o;
     double v;

     va_start(vl, action);
     switch(action) {
     case SMA_NEW: // args: int period
         r.handle = malloc(sizeof(sma_obj));
         r.handle->sma = 0.0;
         r.handle->period = va_arg(vl, int);
         r.handle->values = malloc(r.handle->period * sizeof(double));
         r.handle->lv = 0;
         r.handle->sum = 0.0;
         break;
     case SMA_FREE: // args: sma_obj *handle
         r.handle = va_arg(vl, sma_obj *);
         free(r.handle->values);
         free(r.handle);
         r.handle = NULL;
         break;
     case SMA_VALUES: // args: sma_obj *handle
         o = va_arg(vl, sma_obj *);
         r.values = o->values;
         break;
     case SMA_MEAN: // args: sma_obj *handle
         o = va_arg(vl, sma_obj *);
         r.sma = o->sma;
         break;
     case SMA_ADD: // args: sma_obj *handle, double value
         o = va_arg(vl, sma_obj *);
         v = va_arg(vl, double);
         if ( o->lv < o->period ) {
             o->values[o->lv++] = v;
             o->sum += v;
             o->sma = o->sum / o->lv;
         } else {
             o->sum -= o->values[ o->lv % o->period];
             o->sum += v;
             o->sma = o->sum / o->period;
             o->values[ o->lv % o->period ] = v; o->lv++;
         }
         r.sma = o->sma;
         break;
     }
     va_end(vl);

     return r;
}
