#bin/CuttingStock -b config/testsStageApproach/StdColGenMIP.cfg -a config/dcspParameters.cfg -i data/inputList >> output/StdColGenMIP.out
#mv statistics output/stats_StdColGenMIP
bin/CuttingStock -b config/testsStageApproach/StdColGenCustom.cfg -a config/dcspParameters.cfg -i data/inputList >> output/StdColGenCustom.out
mv statistics output/stats_StdColGenCustom
bin/CuttingStock -b config/testsStageApproach/StgColGenCustom.cfg -a config/dcspParameters.cfg -i data/inputList >> output/StgColGenCustom.out
mv statistics output/stats_StgColGenCustom
