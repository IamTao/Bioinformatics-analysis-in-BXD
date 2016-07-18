HOST=hpc-16
MODE=new

start_cloud:
	gcloud compute instances start ${HOST}

cp_all_cloud:
	gcloud compute copy-files ../../../gene_expression/${MODE}/* ${HOST}:/datadisk/gene_expression
	gcloud compute ssh ${HOST}

cp_code_cloud:
	gcloud compute copy-files ../code/ ${HOST}:/datadisk/gene_expression/
	# gcloud compute ssh ${HOST}

build_platform:
	cd /datadisk/gene_expression/code/
	sh ./.build_visualization.sh
	sh ./.build_computation.sh

run:
	sh .run.sh

backup:
	gcloud compute copy-files ${HOST}:/datadisk/gene_expression/data/updated_genotype/output/ data/updated_genotype/output/

close_cloud:
	gcloud compute instances stop ${HOST}
