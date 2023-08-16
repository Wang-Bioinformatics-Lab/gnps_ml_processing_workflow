run:
	nextflow run ./nf_workflow.nf --resume 

run_docker:
	docker-compose -f docker-compose.yml build
	nextflow run ./nf_workflow.nf --resume -with-docker gnps_ml_processing_workflow-gnps_ml_processing

build:
	docker-compose -D -f docker-compose.yml build

run_interactive:
	docker-compose run --rm gnps_ml_processing bash
