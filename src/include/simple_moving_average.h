#ifndef SMA_H
#define SMA_H

enum Action { SMA_NEW, SMA_FREE, SMA_VALUES, SMA_ADD, SMA_MEAN };

typedef struct sma_obj {
     double  sma;
     double  sum;
     int     period;
     double *values;
     int     lv;
} sma_obj;

typedef union sma_result {
     sma_obj   *handle;
     double     sma;
     double    *values;
} sma_result;

sma_result sma(enum Action action, ...);


#endif /* SMA_H */
