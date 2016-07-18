docker rm -f gene_expression
docker run -it --name gene_expression -v /datadisk/gene_expression/:/home/gene_expression -d itamtao/rstudio /bin/zsh
docker exec -it gene_expression /bin/bash
