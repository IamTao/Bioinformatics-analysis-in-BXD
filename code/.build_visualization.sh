#sudo docker rm -f `sudo docker ps -aq`
docker rm -f shiny
docker run -d --name shiny -p 3838:3838 \
							 -v /datadisk/gene_expression/code/visualization/server:/srv/shiny-server \
						   -v /datadisk/gene_expression/code/visualization/log:/var/log \
						   -v /datadisk/gene_expression/code/visualization/figure:/datadisk/figure \
						   -v /datadisk/gene_expression/data/output:/datadisk/data \
						   --restart=always itamtao/shiny
