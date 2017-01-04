#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "KalmanMuonCalibrator.h"



int
main (int argc, char **argv)
{

  if (argc<=1) {
    printf("No input files specified");
    return 0; 
  }
  KalmanMuonCalibrator *calibrator = new KalmanMuonCalibrator();

  for (int i=1;i<argc;++i) {
    calibrator->processFile(argv[i],"data_13");
  }

  delete calibrator;
  return 0;
}
