HOST=shiny
MODE=new
DISK=datadisk-hpc

start_cloud:
	gcloud compute instances start ${HOST}

cp_all_cloud:
	gcloud compute copy-files ../../../gene_expression/${MODE}/* ${HOST}:/${DISK}/gene_expression
	gcloud compute ssh ${HOST}

cp_code_cloud:
	gcloud compute copy-files ../code ${HOST}:/${DISK}/gene_expression/
	# gcloud compute ssh ${HOST}

build_computation:
	cd /${DISK}/gene_expression/code/
	# build computation
	docker rm -f gene_expression
	docker run -it --name gene_expression -v /${DISK}/gene_expression/:/home/gene_expression -d itamtao/rstudio /bin/zsh
	docker exec -it gene_expression /bin/bash

build_visualization:
	# build shiny
	docker rm -f shiny
	docker run -d --name shiny -p 3838:3838 \
								 -v /${DISK}/gene_expression/code/visualization/server:/srv/shiny-server \
							   -v /${DISK}/gene_expression/code/visualization/log:/var/log \
							   -v /${DISK}/gene_expression/code/visualization/figure:/datadisk/figure \
							   -v /${DISK}/gene_expression/data/output:/datadisk/data \
							   --restart=always itamtao/shiny

run:
	sh .run.sh

backup:
	gcloud compute copy-files ${HOST}:/${DISK}/gene_expression/data/updated_genotype/output/ data/updated_genotype/output/

close_cloud:
	gcloud compute instances stop ${HOST}
