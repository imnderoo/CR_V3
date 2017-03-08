#!/usr/bin/python
# Scripts for converting .BED les to .interval_list les
# BED les are used by the coverage scripts while .interval_list are used by auto-classication

import argparse
import jaydebeapi

def main():
	parser = argparse.ArgumentParser(description='Python version of the legacy createReport.cygwin.sh')
	
	#parser.add_argument('python_module_folder', help = 'Path to metrics.db for storing instrument QC metrics')
	
	args = parser.parse_args()
	
	#lims_db = "/media/sf_PLMGenetcs/Databases/Molecular/Molecular\ Pathology\ Database.accdb"
	lims_db = '/media/sf_PLMGenetics/Development/Molecular/Andrew_Dev/Inherited_Panel/Molecular\\ Pathology\\ Database.accdb'
	#lims_user = "admin"
	#lim_password = "molecularlab"
	
	conn = jaydebeapi.connect('net.ucanaccess.jdbc.UcanaccessDriver', 'jdbc:ucanaccess:///media/sf_PLMGenetics/Development/Molecular/Andrew_Dev/Inherited_Panel/Molecular Pathology Database.accdb;memory=false')

	cursor = conn.cursor()

	query = "SELECT * FROM NGS_QC_Info"
	cursor.execute(query)
	rows = cursor.fetchall()
	for row in rows:
		print row
	cursor.close()
	conn.close()
	
main()
