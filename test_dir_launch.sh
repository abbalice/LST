#####################################
#!/bin/bash

source /soft/centos7/anaconda/2021.05/etc/profile.d/conda.sh
conda activate my_environment

event_loc=Kurils
inpjson=$(cat $event_loc/outputs_step1/set_up.json)
ncellX=$(echo $inpjson | python3 -c "import sys, json; obj=json.load(sys.stdin);print(obj['ncellX'])")
ncellY=$(echo $inpjson | python3 -c "import sys, json; obj=json.load(sys.stdin);print(obj['ncellY'])")

echo $ncellX $ncellY

sed -i s/EVENT_LOC/$event_loc/ test_launch.sh
#sed -i s/NCELLY/$ncellY/ test_launch.sh

start_index=0
end_index=99

if [ $end_index -gt $ncellX  ]; then
  #echo $start_index $end_index
  sed -i s/$start_index-$end_index/$start_index-$ncellX/ test_launch.sh
  qsub -W block=true test_launch.sh  &
else
  qsub -W block=true test_launch.sh  &
fi

while [ $end_index -lt $ncellX ]
do
  #echo $start_index $end_index
  sed -i s/$start_index-$end_index/$(( $start_index + 100 ))-$(( $end_index + 100 ))/ test_launch.sh
  start_index=$(( $start_index + 100 ))
  end_index=$(( $end_index + 100 ))
  if [ $end_index -lt $ncellX  ]; then
     qsub -W block=true test_launch.sh  &
  fi
done

if [ $end_index -gt $ncellX  ]; then
  echo $start_index $end_index
  sed -i s/$start_index-$end_index/$start_index-$ncellX/ test_launch.sh
  qsub -W block=true test_launch.sh  &
fi


sed -i s/$start_index-$ncellX/0-99/ test_launch.sh
sed -i s/$event_loc/EVENT_LOC/ test_launch.sh
#sed -i s/$ncellY/NCELLY/ test_launch.sh

wait

