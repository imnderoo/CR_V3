#/bin/bash

# docker run --rm -v $PWD:/home/ nderoo324/python python /home/GATK_docker_pipeline_160926_gsc_aw.py
# Can't chain docker python  scripts that call more docker scripts.

#AVFKY
#python ./docker_pipeline_20170111b.py /media/sf_PLMGenetics/NextGenData/Validation/CR_V3/170210_M03448_0112_000000000-AVFKY/Data/Intensities/BaseCalls/ /media/sf_resources/createReport_DBs/analysis_type/msh21_inherited/msh21_inherited.bed /media/sf_PLMGenetics/NextGenData/Validation/CR_V3/170210_M03448_0112_000000000-AVFKY/Data/Intensities/BaseCalls/Alignment

python ./create_report.py /media/sf_PLMGenetics/NextGenData/Validation/CR_V3/170210_M03448_0112_000000000-AVFKY/ 20 20 /media/sf_PLMGenetics/NextGenData/Validation/CR_V3/NGS_Analysis/ /media/sf_resources/
