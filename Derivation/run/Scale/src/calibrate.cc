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


  opterr = 0;
  bool MC=false;
  int c;

  while ((c = getopt (argc, argv, "ml:")) != -1)
    switch (c) {
    case 'm':
      MC=true;
      break;
    case 'l':
      calibrator->load(optarg);
      break;
    }



  for (int i=optind;i<argc;++i) {
    if (MC)
      calibrator->processFileMC(argv[i],"tree");
    else
       calibrator->processFile(argv[i],"tree");
    char* newFile;
    if(asprintf(&newFile,"preview_%s",argv[i])<0)
      continue;
    calibrator->save(newFile);
  }

  calibrator->save("calibration.root");
  delete calibrator;
  return 0;
}
