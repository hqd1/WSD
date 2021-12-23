import csv, subprocess
parameter_file_full_path = "/storage/home/hqd1/NeuroProject/parameter.csv"

with open(parameter_file_full_path,"rb") as csvfile:
        reader = csv.reader(csvfile)
        for job in reader:
		print job
            	qsub_command = """qsub -v it={0},id={1},ft={2}, simulation.PBS""".format(*job)
            	exit_status=subprocess.call(qsub_command,shell=True)
            	if exit_status is 1:
                    print "Job {0} failed to submit.".format(qsub_command)
print "Done submitting jobs!"
