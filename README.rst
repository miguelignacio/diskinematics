

H1 Event Shape Analysis (MPIM, 2020)
====================================

.. note::
   Please keep this README file up-to-date

Package structure
-----------
 
  o `Analysis` contains analysis code 

  o `Analysis/Steering` contains steering files for the analysis

  o `LumiFiles` contains the lumi files

  o ... add more here


How to run
----------
A simple start is obtained by:
 
  o Set up H1 software release 4.1.1:  `cd $H1DIST/lcg_releases/releases/4.1.1; . thish1.sh`. 
    Then go back to your work directory

  o Check out repository `git clone https://<user>@stash.desy.de/scm/h1ana/eventshapesmpim.git`

  o Go to analysis directory `cd eventshapesmpim`

  o Compile analysis code `make all`

  o Run a first example with
  ```
  $ ../bin/x86_64-centos7-gcc9-opt/EventShapes -f Steering/es.ep0304dis.steer -n 10000 -o test.root -c DjangoEplus0304_1
  ```

