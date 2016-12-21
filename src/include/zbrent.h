#ifndef ZBRENT_H
#define ZBRENT_H

#include "gday.h"

double zbrent(double (*func) (double, double, double, double, double, double),
              double, double, double, double, double, double, double, double);

#endif /* ZBRENT_H */
