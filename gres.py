#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
import traceback
import os
import argparse
import datetime
import sys
import subprocess
import re
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import shutil
import time
version="1.6.1"
print('gres v{}'.format(version))
start_time = datetime.datetime.now()
call=os.path.abspath(os.getcwd())
contig_list=dict()
element_files=dict()
ERROR=dict()
element_fastas=dict()
old=''
new=''
delimiter=''
contig_size=dict()
not_created=[]
created=[]
genome_ids=dict()
cid_list=[]
contig_id_by_size=dict()
hits_by_qseqid = dict()
organism_names=dict()

parser = argparse.ArgumentParser(
  prog='gres',
  description=
  'Genome Retriever and Element Subtractor',
  epilog='https://github.com/GiulianaPola/gres',add_help=False)
parser.add_argument('-i',help="contig_list file or folder")
parser.add_argument('-o',help='Name of the output folder',default="output_dir")
parser.add_argument('-h', '-help', action='store_true')
parser.add_argument('-v', action='store_true',help='Version')
parser.add_argument('-old',help="Old table with 'Query ID' and 'Organism name' columns")
parser.add_argument('-new',help="New table with 'Organism ID' and 'Organism name' columns")
parser.add_argument('-e',help="Folder with the element's fasta files")
#parser.add_argument('-config', help='Configuration file')
args=parser.parse_args()
def get_actual_path(path,call):
    if " " in path:
      path=path.replace(" ","\ ")
    if path.startswith(".../"):
        actual_path=call.split("/")[:-1]
        actual_path.extend(path.split("/")[1:])
        return os.path.join(actual_path)
    elif os.path.split(path)[0]=='':
        actual_path=call.split("/")
        actual_path.extend(path.split("/"))
        return os.path.join(actual_path)
    else:
        return os.path.abspath(path)
def copy_fasta(exist_fasta,fasta_file,label,contig_id):
    try:
      shutil.copy2(str(exist_fasta), fasta_file)
    except Exception as e:
      log.write("ERROR on copying {} file: {}\n".format(label,str(e)))
      print("ERROR on copying {} file: {}".format(label,str(e)))
      exc_type, exc_value, exc_traceback = sys.exc_info()
      formatted_lines = traceback.format_exc().splitlines()
      for i, line in enumerate(formatted_lines):
        print("{}: {}".format(i, line))
    else:
      log.write("'{}' Copied {} file: {}\n".format(contig_id,label,fasta_file))
def printing_help():
    print('(c) 2023. Giuliana Pola & Arthur Gruber\n')
    print('For more information access: ')
    print('Usage: updateID.py -i <path|file> -old <old>')
    print('Genome Retriever and Element Subtractor')
    print('\nOptional arguments:')
    print('  -h, --help             Show this help message and exit')
    print('  -i <file>        contig_list file with families (required)')
    print('  -o string       Name of the output folder (default: "fasta_files")')
    print('  -old file         Table with Query ID and Organism name columns (required)')
    print('  -new file         Table with Organism ID and Organism name columns (required)')
    print('  -v                    Show version information')
def rename(name):
  i=0
  path = os.path.dirname(name)
  name = os.path.basename(name)
  newname = os.path.join(path, name)
  while os.path.exists(newname):
      i += 1
      newname = os.path.join(path, "{}_{}".format(name, i))
  return newname
def find_delimiter(line):
    delimiters = ['\t', ';',',']
    for delimiter in delimiters:
        if re.match('.*[{}]+.*'.format(delimiter), line):
            return delimiter
    return None
def validate(args):
    args.i=get_actual_path(args.i,call)
    try:
        if os.path.isfile(args.i):
            if os.path.getsize(args.i) > 0:
                with open(args.i, 'r') as file:
                    lines = file.readlines()
                    if len(lines) == 1:
                        delimiter=find_delimiter(lines[0])
                        contig_list[""]=lines.split(delimiter)
                    else:
                        family_index=""
                        contig_index=""
                        delimiter=find_delimiter(lines[0])
                        if delimiter==None:
                          contig_list[""]=lines
                        else:
                          columns=lines[0].replace("\n","").split(delimiter)
                          
                          if "family" in lines[0].lower() or "clade" in lines[0].lower() or "virinae" in lines[0].lower():
                            for i in range(len(columns)):
                              if "family" in columns[i].lower() or "clade" in columns[i].lower() or "virinae" in columns[i].lower():
                                family_index=i
                          if "contig_id" in lines[0].lower() or "contig id" in lines[0].lower() or "query id" in lines[0].lower() or "query_id" in lines[0].lower():
                            for i in range(len(columns)):
                              if "contig_id" in columns[i].lower() or "contig id" in columns[i].lower() or "query id" in columns[i].lower() or "query_id" in columns[i].lower():
                                contig_index=i
                          if family_index=="" and contig_index=="":
                            contig_index=0
                          if len(columns)==1:
                            contig_index=0
                          if contig_index=="" and not family_index=="":
                            if family_index==0:
                              contig_index=1
                            elif family_index==1:
                              contig_index=0
                        for line in lines:
                          if "contig_id" in line.lower() or "contig id" in line.lower() or "query_id" in line.lower() or "query id" in line.lower():
                            pass
                          else:
                            columns=line.replace("\n","").split(delimiter)
                            if not family_index=="":
                              if columns[family_index].isdigit():
                                if not "Family_{}".format(columns[family_index]) in contig_list.keys():
                                  contig_list["Family_{}".format(columns[family_index])]=[]
                                contig_list["Family_{}".format(columns[family_index])].append(columns[contig_index])
                              else:
                                if not columns[family_index] in contig_list.keys():
                                  contig_list[columns[family_index]]=[]
                                contig_list[columns[family_index]].append(columns[contig_index])
                            else:
                              if "Not found" not in contig_list:
                                  contig_list["Not found"] = []  # Initialize as a list if it doesn't exist
                      
                              # Debugging print statements
                              #print("Type of contig_list:", type(contig_list))
                              #print("Type of contig_list['Not found']:", type(contig_list.get("Not found", None)))
                              try:
                                  contig_index = int(contig_index)
                              except ValueError:
                                  pass
                                  #print("Error: contig_index '{}' is not a valid integer.".format(str(contig_index)))
                                  # Handle the error (e.g., set a default value or exit)
                              else:
                                  if contig_index < len(columns):  # Check if contig_index is valid
                                      contig_list["Not found"].append(columns[contig_index])
            else:
                print("ERROR: Contig_id's list (-i) is empty.")
                exit()
        else:
            print("ERROR: Contig_id's list (-i) does not exist.")
            exit()
    except Exception as e:
        print("ERROR: Contig_id's list (-i): ", str(e))
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
           print("{}: {}".format(i, line))
        exit()
    element_folders = []
    if "/*/" in str(args.e):
        path = str(args.e).split("/*/")[0]
        temp_list = str(args.e).split("/*/")[1].split("/")
        while "" in temp_list:
          temp_list.remove("")
        suffix=temp_list[-1]
        for root, dirs, files in os.walk(path):
            for dir in dirs:
                if dir==suffix:
                    element_folders.append(os.path.join(root, dir))
    else:
        element_folders.append(args.e)
    for folder in element_folders:
        folder = get_actual_path(folder, call)
        if os.path.isdir(folder):
            if not folder.endswith("/"):
                folder += "/"
            if not os.path.exists(folder):
                print("ERROR: Element's fasta files' folder '{}' does not exist!".format(folder))
            elif not os.listdir(folder):
                print("ERROR: Element's fasta files' folder '{}' is empty!".format(folder))
            else:
                for root, dirs, files in os.walk(folder):
                    for file in files:
                        if file.endswith('.fasta'):
                            file_path = os.path.join(root, file)
                            contig_id = next((item for v in contig_list.values() for item in v if item in os.path.split(file_path)[-1]),None)
                            if not contig_id==None and contig_id not in element_files:
                                element_files[contig_id] = file_path
                if not element_files:
                    print("ERROR: Input folder '{}' does not contain element's fasta files (.fasta)!".format(folder))
        elif os.path.isfile(folder):
            if not os.path.exists(folder):
                print("ERROR: Element's fasta file '{}' does not exist!".format(folder))
            else:
              element_files[os.path.split(folder)[1]] = folder
    if element_files=={}:
      print("ERROR: No element's fasta file was found in the path (-e)!")
      exit()
    old_lines = []
    args.old=get_actual_path(args.old,call)
    if not os.path.exists(args.old):
        print('ERROR: Old table (-old) does not exist!')
        exit()
    try:
        with open(args.old, mode='r',encoding='utf-8-sig') as f:
            old_lines = f.readlines()
    except IOError:
        print('ERROR: Old table (-old) ccontig_id be read!')
        exit()
    header = old_lines[0].lower()
    delimiter=find_delimiter(header)
    header = header.strip().split(delimiter)
    invalid=['Organism name','Contig ID','Contig size']
    valid=[]
    #header=["Organism ID","Organism name","ContigID_elem#","Contig size","Element coord. on contig","Element size","Gene","Contig coord.","Element coord.","Orientation","Distance(bases)","Gene","Contig coord.","Element coord.","Orientation","Distance(bases)","Gene","Contig coord.","Element coord.","Orientation"]
    for col in ['Organism name','Contig ID','Contig size']:
      found = False
      for c in header:
         if all(substring.lower() in c.lower() for substring in col.split()):
             found = True
             break
      if found:
         invalid.remove(col)
         valid.append(col)
      else:
         print("ERROR: Old table (-old) does not have the '{}' column!".format(col))
    if invalid:
      for col in list(set(invalid)):
           print("ERROR: Old table (-old) does not have the '{}' column!".format(col))
      exit()
    else:
      old=[]
      for line in old_lines:
        old.append(line.strip().split(delimiter))
    new_lines = []
    args.new=get_actual_path(args.new,call)
    if not os.path.exists(args.new):
        print('ERROR: New table (-new) does not exist!')
        exit()
    try:
        with open(args.new, mode='r',encoding='utf-8-sig') as f:
            new_lines = f.readlines()
    except IOError:
        print('ERROR: New table (-new) ccontig_id be read!')
        exit()
    header = new_lines[0].lower()
    delimiter=find_delimiter(header)
    header = header.strip().split(delimiter)
    valid=[]
    invalid=['Genome name','Genome ID']
    for col in ['Genome name','Genome ID']:
      for c in header:
        if all(substring.lower() in c.lower() for substring in col.split()):
            invalid.remove(col)
            valid.append(col)
            break
    if not invalid==[]:
      for col in list(set(invalid)):
          print("ERROR: New table (-new) does not have the '{}' column!".format(col))
      exit()
    else:
      new=[]
      for line in new_lines:
        new.append(line.strip().split(delimiter))
    try:
      if os.path.exists(args.o):
        args.o=rename(args.o)
      os.makedirs(args.o)
      args.o=os.path.abspath(args.o)
    except Exception as e:
      print("ERROR validation:"+str(e))
      print("ERROR: Output folder was not written!")
      exc_type, exc_value, exc_traceback = sys.exc_info()
      formatted_lines = traceback.format_exc().splitlines()
      for i, line in enumerate(formatted_lines):
          print("{}: {}".format(i, line))
      exit()
    return args.i,old,new
def find_file(root_dir, filename):
    for root, dirs, files in os.walk(root_dir):
        if filename in files:
            return os.path.join(root, filename)
    return None
def get_value(table, row,wanted_col, search_col):
    #print("DATETIME: get_value: {}\n".format(datetime.datetime.now()))
    found=False
    columns = table[0]
    search=None
    want=None
    for c in columns:
      if all(sub.lower() in c.lower().split(".")[-1] for sub in wanted_col.split()):
        want = columns.index(c)
      if all(sub.lower() in c.lower().split(".")[-1] for sub in search_col.split()):
        search = columns.index(c)
    if not search==None and not want==None:
      for line in table:
          try:
            if line[search].lower()==row.lower():
              found=True
              return line[want]
              break
          except Exception as e:
            print("'{}' ERROR: get_value exception: {}".format(contig_id,e))
            log.write("'{}' ERROR: get_value exception: {}\n".format(contig_id,e))
            if "ERROR: get_value exception" not in ERROR.keys():
              ERROR['ERROR get_value exception']=[]
            ERROR['ERROR get_value exception'].append(contig_id)              
            exc_type, exc_value, exc_traceback = sys.exc_info()
            formatted_lines = traceback.format_exc().splitlines()
            for i, line in enumerate(formatted_lines):
                print("{}: {}".format(i, line))
      if found==False:
        for line in table:
          if all(sub.lower() in line[search].lower().split() for sub in row.split()):
            found=True
            return line[want]
            break
      if found==False:
        for line in table:
          if row.lower() in line[search].lower():
            found=True
            return line[want]
            break
    return None
def download_sequences(organism_id,output_folder,family):
    #print("DATETIME: download_sequences: {}\n".format(datetime.datetime.now()))
    if not organism_id=="Not found":
      if not family==None:
        genome_folder=os.path.join(output_folder,family,"organism_genomes")
      else:
        genome_folder=os.path.join(output_folder,"organism_genomes")
      if not os.path.isdir(genome_folder):
        os.makedirs(genome_folder)
      os.chdir(genome_folder)
    cmd = "p3-genome-fasta {} > {}.fasta".format(organism_id,organism_id)
    #log.write("'{}' Download the fasta file of sequenced organism genome: {}\n".format(organism_id,cmd))
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print("ERROR '{}': Failed to download organism genome fasta due to {}!\n".format(organism_id, str(e.output)))
        log.write("ERROR '{}': Failed to download organism genome fasta due to {}!\n".format(organism_id, str(e.output)))
        return False
    except Exception as e:
        print("ERROR '{}': Failed to download organism genome fasta due to {}!\n".format(organism_id, e))
        log.write("ERROR '{}': Failed to download organism genome fasta due to {}!\n".format(organism_id, e))
        return False
    else:
        #log.write("'{}' Fasta from organism genome was downloaded!\n".format(organism_id))
        return True
def join_fastas(fasta_files, output_fasta):
    if any(fasta is None for fasta in fasta_files.values()):
        log.write("ATTENTION when joining fastas: One or more input files are None!\n")
        print("ATTENTION when joining fastas: One or more input files are None!")
        keys_to_remove = []
        for contig_id, fasta in fasta_files.items():
            if fasta is None:
                keys_to_remove.append(contig_id)
                log.write("'{}' ATTENTION: Fasta file is None!\n".format(contig_id))
        for key in keys_to_remove:
            fasta_files.pop(key)

    if not all(os.path.isfile(fasta) for fasta in fasta_files.values()):
        log.write("ATTENTION when joining fastas: One or more input files do not exist!\n")
        print("ATTENTION when joining fastas: One or more input files do not exist!")
        keys_to_remove = []
        for contig_id, fasta in fasta_files.items():
            if not os.path.isfile(fasta):
                keys_to_remove.append(contig_id)
                log.write("'{}' ATTENTION: Fasta file do not exist!\n".format(contig_id))
        for key in keys_to_remove:
            fasta_files.pop(key)
    try:
        with open(output_fasta, 'w') as f:
            pass
    except IOError as e:
        log.write("ERROR: Cannot create output file {}!\n".format(output_fasta))
        print("ERROR: Cannot create output file {}!".format(output_fasta))
        log.close()
        exit()
        return False
    if not os.path.exists(output_fasta):
        log.write("ERROR: File {} was not created!\n".format(os.path.basename(output_fasta)))
        print("ERROR: File {} was not created!".format(os.path.basename(output_fasta)))
        log.close()
        exit()
        return False
    try:
        with open(output_fasta, 'w') as output_file:
            for key in fasta_files.keys():
                ffile=fasta_files[key]
                with open(ffile, 'r') as input_file:
                    sequences = []
                    text=input_file.read()
                    text="\n"+text
                    text=text.replace("\n>","\nAQUI COMECA O FASTA>")
                    for fasta in text.split('AQUI COMECA O FASTA'):
                      for line in fasta.split("\n"):
                          if line=="" or line.strip()=="":
                              log.write("")
                          elif line.startswith('>'):
                              line=line.replace(">","")
                              header = [key]
                              header.extend(line.strip().split())
                              header=list(dict.fromkeys(header))
                              header.remove(key)
                              coordinates=""
                              contig=""
                              for value in header:
                                if "-" in value:
                                  coordinates=value
                                  header.remove(coordinates)
                                elif ".con." in value:
                                  key=value
                                  header.remove(value)
                                elif ".con." in text:
                                  if any(c.isalpha() for c in value) and any(c.isdigit() for c in value):
                                    if "element_" not in value.lower():
                                      key=key+".con."+value
                                      header.remove(value)
                              modified_header = '>' + str(key)+" "+coordinates+" "+" ".join(header)
                              sequences.append(modified_header)
                          else:
                              sequences.append(line)
                      for sequence in sequences:
                          output_file.write(sequence)
                          output_file.write("\n")
                    output_file.write("\n")
    except IOError as e:
        log.write("ERROR: File {} was not written!\n".format(os.path.basename(output_fasta)))
        print("ERROR: File {} was not written!".format(os.path.basename(output_fasta)))
    except Exception as e:
        print("ERROR: exception when joining fastas: {}".format(e))
        log.write("ERROR: exception when joining fastas: {}\n".format(e))
        if "ERROR: exception when joining fastas" not in ERROR:
            ERROR["ERROR: exception when joining fastas"] = []
        ERROR["ERROR: exception when joining fastas"].append(os.path.basename(output_fasta))
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
            print("{}: {}".format(i, line))
        log.close()
        exit()
        return False
    if os.path.getsize(output_fasta) == 0:
        log.write("ERROR: File {} is empty!\n".format(os.path.basename(output_fasta)))
        print("ERROR: File {} is empty!".format(os.path.basename(output_fasta)))
        log.close()
        exit()
        return False
    return True
def makeblastdb(query_file,db_file):
  #print("DATETIME: makeblastdb: {}\n".format(datetime.datetime.now()))
  try:
    cmd = "makeblastdb -in {} -out {} -dbtype nucl".format(query_file,db_file)
    log.write("Format the organism genome database: {}\n".format(cmd))
    subprocess.run(cmd, shell=True)
  except Exception as e:
    print("ERROR: exception when joining fastas: {}".format(e))
    log.write("ERROR: exception when joining fastas: {}\n".format(e))
    if "ERROR: exception when joining fastas" not in ERROR.keys():
      ERROR['ERROR: exception when joining fastas']=[]
    ERROR['ERROR: conta_bases.pl exception'].append("")              
    exc_type, exc_value, exc_traceback = sys.exc_info()
    formatted_lines = traceback.format_exc().splitlines()
    for i, line in enumerate(formatted_lines):
        print("{}: {}".format(i, line))
    return False
  else:
    return True
def conta_bases(fasta_file, filtering,contig_id,wanted):
    #print("DATETIME: conta_bases: {}\n".format(datetime.datetime.now()))
    try:
        cmd = "conta_bases.pl -i '{}' |grep '{}'".format(fasta_file, filtering)
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        output, ERROR = process.communicate()
        if ERROR or not output:
            if wanted=="contig_id":
              log.write("'{}' Find contig_id by size:{}\n".format(contig_id,cmd))
            elif wanted=="size":
              log.write("'{}' Find size by contig_id:{}\n".format(contig_id,cmd))
            return None
        elif wanted=="contig_id":
          new_contig_id = output.decode("utf-8").split(">")[-1].split("contig")[0].split("\n")[0].split("\t")[0]
          if "chromosome" in contig_id:
            new_contig_id=contig_id.replace("chromosome","")
          return new_contig_id
        elif wanted=="size":
          size = output.decode("utf-8").split(">")[-1].split("contig")[-1].split("\n")[0].split("\t")[-1]
          return size
    except Exception as e:
        print("'{}' ERROR: conta_bases.pl exception: {}".format(contig_id,e))
        log.write("'{}' ERROR: conta_bases.pl exception: {}\n".format(contig_id,e))
        if "ERROR: conta_bases.pl exception" not in ERROR.keys():
          ERROR['ERROR: conta_bases.pl exception']=[]
        ERROR['ERROR: conta_bases.pl exception'].append(contig_id)              
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
            print("{}: {}".format(i, line))
def find_contig_in_file(new_contig_id, ffile,contig_id):
    #print("DATETIME: find_contig_in_file: {}\n".format(datetime.datetime.now()))
    command = "grep '{}' '{}'".format(new_contig_id, ffile)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    output, err = process.communicate()
    if process.returncode != 0:
        log.write("'{}' Find contig in file: grep '{}' '{}'\n".format(contig_id,new_contig_id, ffile))
        if err is not None:
            log.write("'{}' ERROR: {}\n".format(contig_id, err.decode("utf-8")))
        return False
    elif output:
        return True
        
def count_sequences(filename):
    #print("DATETIME: count_sequences: {}\n".format(datetime.datetime.now()))
    try:
        count = 0
        with open(filename, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                count += 1
        return count
    except Exception as e:
        print("'{}' ERROR: Exception when counting sequences: {}".format(contig_id,e))
        log.write("'{}' ERROR: Exception when counting sequences: {}\n".format(contig_id,e))
        if "ERROR: Exception when counting sequences" not in ERROR.keys():
          ERROR['ERROR: Exception when counting sequences']=[]
        ERROR['ERROR: Exception when counting sequences'].append(contig_id)              
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
            print("{}: {}".format(i, line))
        return None
        
def run_blastn(element_multifasta, organism_multifasta, outfile, outfmt_str):
    #print("DATETIME: run_blastn: {}\n".format(datetime.datetime.now()))
    try:
        blastn_cline = NcbiblastnCommandline(query=element_multifasta, subject=organism_multifasta, out=outfile, outfmt='"7 {}"'.format(outfmt_str), perc_identity=75, num_threads=10)
        log.write("Running BLASTN command: {}\n".format(blastn_cline))
        stdout, stderr = blastn_cline()
    except Exception as e:
        print("ERROR BLASTN: {}".format(e))
        log.write("ERROR BLASTN: {}\n".format(e))
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
            print("{}: {}".format(i, line))
        return False
    else:
        return True
def blastn_elements(element_multifasta, organism_multifasta, outfile):
    organism_id = os.path.basename(organism_multifasta).replace(".fasta","")
    name_columns = { "query id": "qseqid", "query length": "qlen", "q. start": "qstart", "q. end": "qend", "subject id": "sseqid", "subject strand": "sstrand", "subject length": "slen", "s. start": "sstart", "s. end": "send", "% identity": "pident", "alignment length": "length", "mismatches": "mismatch", "gap opens": "gapopen", "evalue": "evalue", "bit score": "bitscore", "subject seq": "sseq" }
    outfmt_str = "qseqid qlen qstart qend sseqid sstrand slen sstart send pident length mismatch gapopen evalue bitscore"
    cols = []
    families = []
    extraction=dict()
    qseqid=""
    summary=dict()
    blast_contigs=[]
    try:
        os.chdir(os.path.join(args.o,"BLASTN"))
        if run_blastn(element_multifasta, organism_multifasta, outfile, outfmt_str):
            if os.path.isfile(outfile):
                if os.path.getsize(outfile) > 0:
                    fields = []
                    with open(outfile, 'r') as f:
                        for line in f.readlines():
                            hit = dict()
                            if "#" in line:
                                if "Query:" in line:
                                  qseqid = line.split("Query: ")[-1].split()[0]
                                  #if qseqid not in summary.keys():
                                    #summary[qseqid] = dict()
                                if "Fields:" in line and fields == []:
                                    fields = line.replace("\n","").replace(": ",":").replace(", ",",").split(":")[-1].split(",")
                                    for field in fields:
                                        cols.append(name_columns[field])
                            else:
                                category=[]
                                pident = line.split()[fields.index("% identity")]
                                lenght = line.split()[fields.index("alignment length")]
                                qlenght = line.split()[fields.index("query length")]
                                perclength=int(lenght)/int(qlenght)*100
        
                                for col, value in zip(cols, line.split()):
                                    hit[col] = value
                                sseqid = line.split()[fields.index("subject id")]
                                qseqid = line.split()[fields.index("query id")]
                                if float(pident) < 75 or int(lenght) < 1000:
                                    if float(pident) < 75:
                                      category.append("pident<{}".format(str(75)))
                                    if float(pident) < 1000:
                                      category.append("alignment lenght<{}".format(str(1000)))
                                    if float(perclength) < 98:
                                      category.append("alignment lenght<{}".format(str(1000)))
                                    category="Invalid hit ("+" and ".join(category)+")"
                                elif float(pident) >= 75 and int(lenght) >= 1000:
                                    if qseqid not in blast_contigs:
                                      blast_contigs.append(qseqid)
                                    family = next((f for f, contigs in contig_list.items() if qseqid in contigs), None)
                                    if not family in families:
                                        families.append(family)
                                    if float(pident) >= 98:
                                        category.append("pident>={}%".format(str(98)))
                                    if float(lenght) >= 3000:
                                        category.append("alignment lenght>={}".format(str(3000)))
                                    if float(perclength) >= 98:
                                        category.append("% alignment lenght>={}%".format(str(98)))
                                      category="Valid hit for subtraction ("+" and ".join(category)+")"
                                    if float(pident) >= 98 and int(lenght) >= 3000 and perclength>=98:
                                      category=[]
                                      if organism_id not in extraction.keys():
                                          extraction[organism_id] = dict()
                                      if family not in extraction[organism_id].keys():
                                          extraction[organism_id][family] = []
                                      if qseqid not in extraction[organism_id][family]:
                                          extraction[organism_id][family].append(qseqid)
                                      if float(pident) >= 98:
                                        category.append("pident>={}%".format(str(98)))
                                      if float(lenght) >= 3000:
                                        category.append("alignment lenght>={}".format(str(3000)))
                                      if float(perclength) >= 98:
                                        category.append("% alignment lenght>={}%".format(str(98)))
                                      category="Valid hit for extraction ("+" and ".join(category)+")"
                                      
                                      if qseqid not in hits_by_qseqid.keys():
                                          hits_by_qseqid[qseqid] = dict()
                                      if organism_id not in hits_by_qseqid[qseqid].keys():
                                          hits_by_qseqid[qseqid][organism_id] = []
                                      if hit not in hits_by_qseqid[qseqid][organism_id]:
                                          hits_by_qseqid[qseqid][organism_id].append(hit)
                                      with open(os.path.join(args.o,"selected_new_names.tab"), 'r+') as f:
                                        lines = f.readlines()
                                        found = False
                                        for i, line in enumerate(lines):
                                            columns = line.replace("\n", "").split("\t")
                                            if "Family" in line and "Old" in line and "New" in line and "Organism" in line and "Method" in line:
                                              lines.remove(line)
                                            if family in line and hit["qseqid"] in line and hit["sseqid"] in line and organism_id in line:
                                                if not "Blast" in line:
                                                    columns[0] = family
                                                    columns[1] = hit["qseqid"]
                                                    columns[2] = hit["sseqid"]
                                                    columns[3] = organism_id
                                                    columns[4] = organism_names[organism_id]
                                                    columns[5] = columns[5].replace(" ","").split(",")
                                                    columns[5].append("Blast")
                                                    columns[5] = ",".join(sorted(list(set(columns[5]))))
                                                    lines[i] = "\t".join(columns) + '\n'
                                                found = True
                                                break
                                        if not found:
                                            lines.append("{}\t{}\t{}\t{}\t{}\tBlast\n".format(family, hit["qseqid"], hit["sseqid"], organism_id,organism_names[organism_id]))
                                        max_lengths = [max(len(column) for column in [line.split("\t")[i] for line in lines]) for i in range(6)]
                                        f.seek(0)
                                        f.write("{:<{}}\t{:<{}}\t{:<{}}\t{:<{}}\t{:<{}}\t{:<{}}\n".format("Family", max_lengths[0], "Old", max_lengths[1], "New", max_lengths[2], "Organism ID", max_lengths[3], "Organism name", max_lengths[4], "Method", max_lengths[5]))
                                        for line in sorted(list(set(lines))):
                                            columns = line.replace("\n", "").split("\t")
                                            f.write("{:<{}}\t{:<{}}\t{:<{}}\t{:<{}}\t{:<{}}\t{:<{}}\n".format(columns[0], max_lengths[0], columns[1], max_lengths[1], columns[2], max_lengths[2], columns[3], max_lengths[3], columns[4], max_lengths[4], columns[5], max_lengths[5]))
                                        f.truncate()
                                    else:
                                      if float(pident) >= 75:
                                        category.append("pident>={}".format(str(75)))
                                      if float(lenght) >= 1000:
                                        category.append("lenght>={}".format(str(1000)))
                                      category="Valid hit for subtraction ("+" and ".join(category)+")"
                                    if sseqid not in hits_by_sseqid.keys():
                                        hits_by_sseqid[sseqid] = []
                                    if hit not in hits_by_sseqid[sseqid]:
                                        hits_by_sseqid[sseqid].append(hit)
                                    if organism_id not in genome_ids.keys():
                                        genome_ids[organism_id] = dict()
                                    if family not in genome_ids[organism_id].keys():
                                        genome_ids[organism_id][family] = []
                                    if qseqid not in genome_ids[organism_id][family]:
                                        genome_ids[organism_id][family].append(qseqid)
                                if category not in summary.keys():
                                    summary[category] = []
                                if qseqid not in summary[category]:
                                    summary[category].append(qseqid)  

                        if organism_id in extraction.keys():
                          if len(extraction[organism_id].keys())>1:
                            print("ATTENTION: Organism's genome ({}) has more than one family ({}) of elements ({})!".format( organism_id, ",".join(extraction[organism_id].keys()), ",".join(str(value) for value in extraction[organism_id].values())))
                            log.write("ATTENTION: Organism's genome ({}) has more than one family ({}) of elements ({})!\n".format( organism_id, ",".join(extraction[organism_id].keys()), ",".join(str(value) for value in extraction[organism_id].values())))
                

#                    for wrong_family in sorted(list(set(contig_list.keys()))):
#                        fasta_file=os.path.join(args.o,wrong_family,"organism_genomes","{}.fasta".format(organism_id))
#                        if not wrong_family in families:
#                            if os.path.isfile(fasta_file):
#                                os.remove(fasta_file)
                    for contig_id in sorted(list(set(blast_contigs))):
                            family = next((f for f, contigs in contig_list.items() if contig_id in contigs), None)
                            fasta_file=os.path.join(args.o,family,"organism_genomes","{}.fasta".format(organism_id))
                            exist_fasta=find_file(args.o, "{}.fasta".format(organism_id))
                            if os.path.isfile(exist_fasta) and not os.path.isfile(fasta_file):
                                copy_fasta(exist_fasta,fasta_file,"organism genome",contig_id)
                    log.write("'{}' BLASTN Summary: {}\n".format(organism_id, str(summary)))
                    if hits_by_qseqid==dict() and hits_by_sseqid==dict():
                       log.write("ATTENTION: Organism genome '{}' has no valid hits (pident>={}% and alignment lenght>={})!\n".format(organism_id,str(75),str(1000)))
                       print("ATTENTION: Organism genome '{}' has no valid hits (pident>={}% and alignment lenght>={})!".format(organism_id,str(75),str(1000)))
                       return hits_by_qseqid,hits_by_sseqid,genome_ids
                    else:
                      return hits_by_qseqid,hits_by_sseqid,genome_ids
                else:
                  log.write("'{}' ERROR BLASTN: No hits found!\n".format(organism_id))
                  print("'{}' ERROR BLASTN: No hits found!".format(organism_id))
                  return hits_by_qseqid,hits_by_sseqid,genome_ids
            else:
                log.write("'{}' ERROR BLASTN: Result file was not created!\n".format(organism_id))
                print("'{}' ERROR BLASTN: Result file was not created!".format(organism_id))
                return hits_by_qseqid,hits_by_sseqid,genome_ids
        else:
              log.write("'{}' ERROR BLASTN: Execution failed!\n".format(organism_id))
              print("'{}' ERROR BLASTN: Execution failed!".format(organism_id))
              return hits_by_qseqid,hits_by_sseqid,genome_ids
    except Exception as e:
          log.write("ERROR BLASTN: {}\n".format(str(e)))
          print("ERROR BLASTN: {}".format(str(e)))
          exc_type, exc_value, exc_traceback = sys.exc_info()
          formatted_lines = traceback.format_exc().splitlines()
          for i, line in enumerate(formatted_lines):
              log.write("{}: {}\n".format(i, line))
          exit()
          return hits_by_qseqid,hits_by_sseqid,genome_ids
def extract_element(organism_id, fasta_folder, new_contig_id,old_contig_id,start,end,sstrand):
    #print("DATETIME: extract_element: {}\n".format(datetime.datetime.now()))
    try:
      found=False
      output_file = os.path.join(fasta_folder, "element", old_contig_id + '.fasta')
      with open(output_file, 'w') as outfile:
          fasta_file = os.path.join(fasta_folder,"organism_genomes", organism_id + '.fasta')
          with open(fasta_file, 'r') as fasta_handle:
              for record in SeqIO.parse(fasta_handle, 'fasta'):
                  if record.id == new_contig_id:
                      found=True
                      if sstrand=="plus":
                          seq = str(record.seq)[int(start)-1:int(end)]
                      elif sstrand=="minus":
                          seq = str(Seq(str(record.seq)[int(end)-1:int(start)]))
                          #.reverse_complement()?
                      new_record = SeqRecord(Seq(seq), id=new_contig_id , description="{}-{}".format(str(start),str(end)))
                      SeqIO.write(new_record, outfile, "fasta")
                      sseq[old_contig_id]=seq
                      break
      if found==True:
        return True
      else:
        print("'{}' ERROR: Contig {} not found in FASTA file!".format(old_contig_id,new_contig_id))
        log.write("'{}' ERROR: Contig {} not found in FASTA file!\n".format(old_contig_id,new_contig_id))
        return False
    except Exception as e: 
        print("'{}' ERROR: {}".format(old_contig_id, str(e)))
        log.write("'{}' ERROR: {}\n".format(old_contig_id, str(e)))
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
          print("{}: {}".format(i, line))
        return False
    else:
        return True          
def write_records_to_fasta(records, output_file):
    #print("DATETIME: write_records_to_fasta: {}\n".format(datetime.datetime.now()))
    try:
        dir_path = os.path.dirname(output_file)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        with open(output_file, "w") as f:
            SeqIO.write(records, f, "fasta")
    except Exception as e: 
        print("'{}' ERROR: {}".format(organism_id, str(e)))
        log.write("'{}' ERROR: {}\n".format(organism_id, str(e)))
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
            print("{}: {}".format(i, line))
        return False
    else:
        return True
def subtract_element(genome_file, output_file, hits_by_sseqid):
    organism_id = os.path.basename(genome_file).replace(".fasta", "")
    try:
        records = list(SeqIO.parse(genome_file, "fasta"))
        changed = False
        sublog=open(os.path.join(args.o,"subtraction.log"), 'a')
        for rec in records:
            if str(rec.id) in list(hits_by_sseqid.keys()):
                
                coordinates=[]
                for hit in hits_by_sseqid[rec.id]:
                  coordinates.append((int(hit["sstart"]), int(hit["send"])))
                  sublog.write(str(hit) + '\n')
                coordinates = list(set(coordinates))
                coordinates.sort(key=lambda x: x[0])
                new_coordinates = []
                for start, end in coordinates:
                    if end < start:
                        start, end = end, start
                    start -= 1
                    end -= 1
                    new_coordinates.append((start, end))        
                coordinates = new_coordinates
                coordinates=sorted(list(set(coordinates)))
                log.write("'{}': '{}' Coordinates: [{}]\n".format(organism_id,rec.id,','.join(map(str, sorted(coordinates)))))
                gaps = set(range(len(rec.seq)))
                for start, end in coordinates:
                    element = set(range(start, end))
                    gaps -= element
                pairs = []
                i = 0
                gaps=list(gaps)
                while i < len(gaps):
                    start = gaps[i]
                    while i + 1 < len(gaps) and gaps[i + 1] == gaps[i] + 1:
                        i += 1
                    if gaps[i] != start:
                        pairs.append((start, gaps[i]))
                    i += 1
                log.write("'{}': '{}' Gaps: [{}]\n".format(organism_id,rec.id,','.join(map(str, pairs))))
                #new_seq = ''.join([str(rec.seq)[nucl] for nucl in gaps])
                if not gaps==[]:
                    changed = True
                    new_seq = "".join(str(rec.seq)[s:e] for s, e in pairs)
                    rec.seq = Seq(new_seq)
        sublog.close()
        if changed:
            if not write_records_to_fasta(records, output_file):
                return False
    except Exception as e:
        log.write("'{}' ERROR: {}\n".format(organism_id, str(e)))
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
            log.write("{}: {}\n".format(i, line))
        return False
    else:
        log.write("'{}': Successfully subtracted genome file: {}\n".format(organism_id,output_file))
        return True
if len(sys.argv)==1:
    printing_help()
    exit()
if "-h" in sys.argv or "--help" in sys.argv:
    printing_help()
    exit()
if args.v:
    print(version)
    exit()
else:
  valid=True
  if not args.i:
    print("ERROR: Missing contig_id's list (-i)!")
    valid=False
  if not args.old:
    print("ERROR: Missing old table with 'Query ID' and 'Organism name' columns (-old)!")
    valid=False
  if not args.new:
    print("ERROR: Missing new table with 'Organism ID' and 'Organism name' columns (-new)!")
    valid=False
  if not args.e:
    print("ERROR: Missing folder with the element's fasta files (-e)!")
    valid=False
  if valid==False:
    exit()
arg_input,old,new = validate(args)
try:
  log=open(os.path.join(args.o,"file.log"), 'w')
  try:
    log.write('gres v{}\n'.format(version))
    log.write('(c) 2023. Giuliana Pola & Arthur Gruber\n')
    log.write('For more information access: ')
    log.write('\nStart time: {}\n'.format(start_time.strftime("%d/%m/%Y, %H:%M:%S")))
    #print("DATETIME: gres.py: {}\n".format(datetime.datetime.now()))
    log.write('\nWorking directory: {}\n'.format(call))
    log.write('\nCommand line: {}\n'.format(' '.join(sys.argv)))
    sublog=open(os.path.join(args.o,"subtraction.log"), 'w')
    sublog.close()
    user=""
    try:
      user=os.getlogin()
    except Exception as e:
      try:
        user=os.environ['LOGNAME']
      except Exception as e:
        try:
          user=os.environ['USER']
        except Exception as e:
          pass
        else:
          pass
      else:
        pass
    else:
      pass
    if not user=="":
      log.write('\nUser: {}\n'.format(user))
    log.write('\nParameters:\n')
    for arg in vars(args):
        value = getattr(args, arg)
        if value is not None and value is not False:
            log.write("{}={}\n".format(arg,value))
  except Exception as e:
    print("ERROR: Log file was not written!")
    exit()
except Exception as e:
    print("ERROR: Log file was not created!")
    exit()
log.close()
try:
  new_names=open(os.path.join(args.o,"selected_new_names.tab"), 'w')
except Exception as e:
  print("ERROR: Selected_new_names.tab file was not created!")
  exit()
else:
  try:
    new_names.write("Family\tOld\tNew\tOrganism ID\tOrganism name\tMethod\n")
  except Exception as e:
    print("ERROR: Selected_new_names.tab file was not written!")
    exit()
new_names.close()
genomes=dict()
genome_ids=dict()
genomes_files=dict()
element_fastas=dict()
with open(os.path.join(args.o, "file.log"), 'a') as log:
    #log.write("\nELEMENT_FILES={}\n".format(element_files))
    for family in sorted(list(set(contig_list.keys()))):
        log.write("\n{}\n".format(family.upper()))
        fasta_folder = os.path.join(args.o, family) if family else args.o
        if not os.path.exists(fasta_folder):
            os.makedirs(fasta_folder.replace("\ ", " "))
        contig_list[family] = list(set(contig_list[family]))
        for contig_id in sorted(contig_list[family]):
            log.write("\nCONTIG_ID={}\n".format(contig_id.upper()))
            element_fastas[contig_id] = next((element_files[k] for k in element_files.keys() if contig_id in k), None)
            if os.path.split(args.old)[-1] == "":
                tail = "table (-old)"
            else:
                tail = os.path.split(args.old)[-1] + " (-old)"
                
            if not find_contig_in_file("_".join(contig_id.split("_")[:-1]), args.old, "_".join(contig_id.split("_")[:-1])):
                log.write("'{}' ERROR: Contig_id '{}' not found in {}!\n".format(contig_id,"_".join(contig_id.split("_")[:-1]), tail))
                print("'{}' ERROR: Contig_id '{}' not found in {}!".format(contig_id,"_".join(contig_id.split("_")[:-1]), tail))
                if "ERROR: Contig_id not found in {}".format(tail) not in ERROR:
                    ERROR["ERROR: Contig_id not found in {}".format(tail)] = []
                ERROR["ERROR: Contig_id not found in {}".format(tail)].append(contig_id)
            else:
                if not find_contig_in_file(contig_id, args.old, contig_id):
                  log.write("'{}' ERROR: Contig_id not found in {}!\n".format(contig_id, tail))
                  print("'{}' ERROR: Contig_id not found in {}!".format(contig_id, tail))
                  if "ERROR: Contig_id not found in {}".format(tail) not in ERROR:
                      ERROR["ERROR: Contig_id not found in {}".format(tail)] = []
                  ERROR["ERROR: Contig_id not found in {}".format(tail)].append(contig_id)
                  organism_name = get_value(old, "_".join(contig_id.split("_")[:-1]), "Organism name", "contig id")
                else:
                  organism_name = get_value(old, contig_id, "Organism name", "contig id")
                if organism_name is None:
                    log.write("ERROR '{}': Organism name was not found in {}!\n".format(contig_id, tail))
                    print("ERROR '{}': Organism name was not found in {}!".format(contig_id, tail))
                    if "ERROR: Organism_name was not found in {}".format(tail) not in ERROR:
                        ERROR["ERROR: organism_name was not found in {}".format(tail)] = []
                    ERROR["ERROR: Organism_name was not found in {}".format(tail)].append(contig_id)
                else:
                    log.write("'{}' Organism_name ({}): {}\n".format(contig_id, tail, organism_name))
                    values_old = {}
                    element_size = get_value(old, "_".join(contig_id.split("_")[:-1]), "Element size", "contig id")
                    contig_size = get_value(old, "_".join(contig_id.split("_")[:-1]), "contig size", "contig id")
                    old_organism_id = get_value(old, "_".join(contig_id.split("_")[:-1]), "Organism id", "contig id")
                    organism_id = None
                    if element_size is not None:
                        values_old["element_size"] = element_size
                    if contig_size is None:
                        log.write("ERROR '{}': contig_size was not found in {}!\n".format(contig_id, tail))
                        print("ERROR '{}': contig_size was not found in {}!".format(contig_id, tail))
                        if "ERROR: contig_size was not found in {}".format(tail) not in ERROR:
                            ERROR["ERROR: contig_size was not found in {}".format(tail)] = []
                        ERROR["ERROR: contig_size was not found in {}".format(tail)].append(contig_id)
                    else:
                        values_old["contig_size"] = contig_size
                        if old_organism_id is not None:
                            log.write("'{}' Organism_id ({}): {}\n".format(contig_id, tail, old_organism_id))
                            values_old["organism_id"] = old_organism_id
                            if old_organism_id not in organism_names.keys():
                              organism_names[old_organism_id]=organism_name
                    if len(values_old) > 0:
                        log.write("'{}' {}: {}\n".format(contig_id, tail, values_old))
                    if organism_name is not None:
                        organism_id = get_value(new, organism_name, "genome id", "genome name")
                        if os.path.split(args.new)[-1] == "":
                            tail = "table -new"
                        else:
                            tail = os.path.split(args.new)[-1] + " (-new)"
                    if organism_id==None:
                        log.write("'{}' find organism_id in file: grep '{}' '{}'\n".format(contig_id, organism_name, args.new))
                        log.write("ERROR '{}': Organism_id was not found in {}!\n".format(contig_id, tail))
                        print("ERROR '{}': Organism_id was not found in {}!".format(contig_id, tail))
                        if not "ERROR: Organism_id was not found in {}".format(tail) in ERROR.keys():
                          ERROR["ERROR: Organism_id was not found in {}".format(tail)]=[]
                        ERROR["ERROR: Organism_id was not found in {}".format(tail)].append(contig_id)
                    else:
                        if os.path.split(args.new)[-1]=="":
                          tail="-new"
                        else:
                          tail=os.path.split(args.new)[-1]+" (-new)"
                        log.write("'{}' New organism_id ({}): {}\n".format(contig_id,tail,organism_id))
                        if organism_id not in organism_names.keys():
                          organism_names[organism_id]=organism_name
                    if organism_id is not None and old_organism_id is not None:
                        if old_organism_id!= organism_id:
                          log.write("'{}' ATTENTION: organism_id has changed!\n".format(contig_id))
                    list_organism_id=[]
                    #if not old_organism_id==None:
                          #list_organism_id.append(old_organism_id)
                    if not organism_id==None:
                          list_organism_id.append(organism_id)
                    else:
                          list_organism_id.append(old_organism_id)
                    
                    for key_name in list_organism_id:
                        if key_name not in genome_ids.keys():
                            genome_ids[key_name]=dict()
                        if family not in genome_ids[key_name]:
                          genome_ids[key_name][family] = []
                        if contig_id not in genome_ids[key_name][family]:
                          genome_ids[key_name][family].append(contig_id)
                        family=next((fam for fam, contigs in contig_list.items() if contig_id in contigs), None)
                        if not family==None:
                          genome_folder=os.path.join(args.o,family,"organism_genomes")
                        else:
                          genome_folder=os.path.join(args.o,"organism_genomes")
                        if not os.path.isdir(genome_folder):
                          os.makedirs(genome_folder)
                        os.chdir(genome_folder)
                        genome_fasta=os.path.join(genome_folder, "{}.fasta".format(key_name))
                        if organism_id in genomes_files.keys():
                          exist_fasta=genomes_files[organism_id]
                        elif old_organism_id in genomes_files.keys():
                          exist_fasta=genomes_files[old_organism_id]
                        else:
                          exist_fasta=[]
                        if family in genome_ids[key_name].keys():
                           if not exist_fasta==[]:
                             if os.path.isfile(exist_fasta) and not os.path.isfile(genome_fasta):
                               copy_fasta(exist_fasta,genome_fasta,"organism genome",contig_id)
                           else:
                               if not download_sequences(key_name,args.o,family):  
                                 if not "ERROR: Sequenced genomes of the organism have not been downloaded" in ERROR.keys():
                                   ERROR["ERROR: Sequenced genomes of the organism have not been downloaded"]=[]
                                 ERROR["ERROR: Sequenced genomes of the organism have not been downloaded"].append(contig_id)
                                 log.write("'{}' ERROR: Sequenced genomes of the organism have not been downloaded!\n".format(contig_id))
                                 print("'{}' ERROR: Sequenced genomes of the organism have not been downloaded!".format(contig_id))
                               else:
                                 log.write("'{}' Organism genome file: {}\n".format(contig_id,genome_fasta))
                           if not os.path.exists(os.path.join(genome_folder,"{}.fasta".format(key_name))):
                                 if not "ERROR: Fasta file of sequenced genomes of the organism was not created" in ERROR.keys():
                                   ERROR["ERROR: Fasta file of sequenced genomes of the organism was not created"]=[]
                                 ERROR["ERROR: Fasta file of sequenced genomes of the organism was not created"].append(contig_id)
                                 log.write("'{}' ERROR: Fasta file of sequenced genomes of the organism was not created!: {}\n".format(contig_id,os.path.join(genome_folder,"{}.fasta".format(key_name))))
                                 print("'{}' ERROR: Fasta file of sequenced genomes of the organism was not created!: {}".format(contig_id,os.path.join(genome_folder,"{}.fasta".format(key_name))))
                           elif not os.path.getsize(os.path.join(genome_folder,"{}.fasta".format(key_name))) > 0:
                                 if not "ERROR: Fasta file of sequenced genomes of the organism is empty" in ERROR.keys():
                                   ERROR["ERROR: Fasta file of sequenced genomes of the organism is empty"]=[]
                                 ERROR["ERROR: Fasta file of sequenced genomes of the organism is empty"].append(contig_id)
                                 log.write("'{}' ERROR: Fasta file of sequenced genomes of the organism is empty!\n".format(contig_id))
                                 print("'{}' ERROR: Fasta file of sequenced genomes of the organism is empty!".format(contig_id))
                           else:
                                 log.close()
                                 log=open(os.path.join(args.o,"file.log"), 'a')
                                 if key_name not in genomes_files.keys():
                                   genomes_files[key_name]=os.path.join(genome_folder,"{}.fasta".format(key_name))
                                 #log.write("'{}' Organism's genome fasta: {}\n".format(contig_id,os.path.join(genome_folder,"{}.fasta".format(key_name))))   
                        if family not in genomes.keys():
                          genomes[family]=dict()
                        if key_name not in genomes[family]:
                          genomes[family][key_name] = []
                        if contig_id not in genomes[family][key_name]:
                          genomes[family][key_name].append(contig_id)
                        genome_file=os.path.join(genome_folder,"{}.fasta".format(key_name))
                        if not os.path.isfile(genome_file):
                          if key_name in genomes_files.keys():
                            genomes_file=os.path.join(genome_folder,"{}.fasta".format(key_name))
                          else:
                            genomes_file=find_file(args.o, "{}.fasta".format(key_name))     
                    log.close()
                    log=open(os.path.join(args.o,"file.log"), 'a')
#log.write("\nELEMENT_FASTAS={}\n".format(element_fastas))
if join_fastas(element_fastas, os.path.join(args.o,"all_element.fasta")):
  for organism_id in list(set(genome_ids.keys())):
    log.close()
    log=open(os.path.join(args.o,"file.log"), 'a')
    sseq=dict()
    hits_by_sseqid=dict()
    log=open(os.path.join(args.o,"file.log"), 'a')
    log.write("\nORGANISM_ID: {}\n".format(organism_id.upper()))
    if not os.path.isdir(os.path.join(args.o,"BLASTN")):
      os.makedirs(os.path.join(args.o,"BLASTN"))
    os.chdir(os.path.join(args.o,"BLASTN"))
    genome_ids=dict()
    element_multifasta=os.path.join(args.o,"all_element.fasta")
    organism_multifasta=find_file(args.o, "{}.fasta".format(organism_id)) 
    outfile=os.path.join(args.o,"BLASTN","blastn_{}.txt".format(organism_id))
    hits_by_qseqid,hits_by_sseqid,genome_ids=blastn_elements(element_multifasta, organism_multifasta,outfile)
    
    #genome_ids[organism_id][family].append(qseqid)
    #log.write("genome_ids='{}'\n".format(genome_ids))
    #log.write("\n'{}' hits={}\n".format(organism_id,hits_by_qseqid))
    log.close()
    log=open(os.path.join(args.o,"file.log"), 'a')
    if hits_by_sseqid==None:
      if not "ERROR: No hits found or valid" in ERROR.keys():
        ERROR["ERROR: No hits found or valid"]=[]
      ERROR["ERROR: No hits found or valid"].append(organism_id)
      print("'{}' ERROR: No hits found or valid".format(organism_id))
      log.write("'{}' ERROR: No hits found or valid!\n".format(organism_id))
    else:
      family=next((fam for fam, contigs in contig_list.items() if contig_id in contigs), None)
      if not family==None:
        genome_folder=os.path.join(args.o,family,"organism_genomes")
      else:
        genome_folder=os.path.join(args.o,"organism_genomes")
      fasta_file = genomes_files[organism_id]
      if not os.path.isfile(fasta_file):
        download_sequences(organism_id,genome_folder,family)
      nseqs=count_sequences(fasta_file)
      os.chdir(genome_folder)
      if nseqs==None:
        if not "ERROR: Couldn't count the contigs in the organism's multifasta" in ERROR.keys():
          ERROR["ERROR: Couldn't count the contigs in the organism's multifasta"]=[]
        ERROR["ERROR: Couldn't count the contigs in the organism's multifasta"].extend(genome_ids[organism_id].values())
        log.write("'{}' ERROR: Couldn't count the contigs in the Fasta file of the organism's multifasta!\n".format(organism_id))
        print("'{}' ERROR: Couldn't count the contigs in the Fasta file of the organism's multifasta!".format(organism_id))
      elif nseqs==0:
        if not "ERROR: The organism's multifasta has no contig" in ERROR.keys():
          ERROR["ERROR: The organism's multifasta has no contig"]=[]
        ERROR["ERROR: The organism's multifasta has no contig"].extend(genome_ids[organism_id].values())
        log.write("'{}' ERROR: The organism's multifasta has no contig!\n".format(organism_id))
        print("'{}' ERROR: Couldn't count the contigs in the Fasta file of the organism's multifasta!".format(organism_id))
        for family in list(set(genome_ids[organism_id].keys())):
          for contig_id in list(set(genome_ids[organism_id][family])):
            if contig_id in created:
              created.remove(contig_id)
            if contig_id not in not_created:
              not_created.append(contig_id)
      else:
        if nseqs==1:
          if not "ATTENTION: The organism's multifasta has only one contig" in ERROR.keys():
            ERROR["ATTENTION: The organism's multifasta has only one contig"]=[]
          ERROR["ATTENTION: The organism's multifasta has only one contig"].append(organism_id)
          log.write("'{}' ATTENTION: The organism's multifasta has only one contig!\n".format(organism_id))
        if not organism_id in genome_ids.keys():
           log.write("Organism genome '{}' has no valid hits!\n".format(organism_id))
        else:
          family=sorted(list(set(genome_ids[organism_id].keys())))[0]
          if not family=="":
            fasta_folder=os.path.join(args.o,family)
          else:
            fasta_folder=args.o
          #subtracted_file=os.path.join(fasta_folder, "without_element", "{}-{}.fasta".format(organism_id,",".join(hits_by_qseqid.keys())))
          subtracted_file=os.path.join(fasta_folder, "without_element", "{}_without_element.fasta".format(organism_id))
          genomes_file=find_file(args.o,"{}.fasta".format(organism_id))
          if not subtract_element(genomes_file, subtracted_file, hits_by_sseqid):
            if not "ERROR: Exception when subtracting contig" in ERROR.keys():
              ERROR["ERROR: Exception when subtracting contig"]=[]
            ERROR["ERROR: Exception when subtracting contig"].append(organism_id)
            log.write("'{}' ERROR: Exception when subtracting contig!\n".format(organism_id))
            print("'{}' ERROR: Exception when subtracting contig!".format(organism_id))
            for family in list(set(genome_ids[organism_id].keys())):
              for contig_id in list(set(genome_ids[organism_id][family])):
                if contig_id not in not_created:
                    not_created.append(contig_id)
                if contig_id in created:
                    created.remove(contig_id)
          elif not os.path.exists(subtracted_file):
            log.write("'{}' ERROR: Fasta file of the subtracted contig was not created!\n".format(organism_id,",".join(genome_ids[organism_id])))
            print("'{}' ERROR: Fasta file of the subtracted contig was not created!".format(organism_id,",".join(genome_ids[organism_id]))) 
            if not "ERROR: Fasta file of the subtracted contig was not created" in ERROR.keys():
              ERROR["ERROR: Fasta file of the subtracted contig was not created"]=[]
            ERROR["ERROR: Fasta file of the subtracted contig was not created"].extend(genome_ids[organism_id])
            for family in list(set(genome_ids[organism_id].keys())):
              for contig_id in list(set(genome_ids[organism_id][family])):
                if contig_id not in not_created:
                    not_created.append(contig_id)
                if contig_id in created:
                    created.remove(contig_id)
          elif not os.path.getsize(subtracted_file) > 0:
            log.write("'{}' ERROR: Fasta file of the subtracted contig is empty!\n".format(organism_id,",".join(genome_ids[organism_id])))
            print("'{}' ERROR: Fasta file of the subtracted contig is empty!".format(organism_id,",".join(genome_ids[organism_id])))  
            if not "ERROR: Fasta file of the subtracted contig is empty" in ERROR.keys():
              ERROR["ERROR: Fasta file of the subtracted contig is empty"]=[]
            ERROR["ERROR: Fasta file of the subtracted contig is empty"].extend(genome_ids[organism_id])
            for family in list(set(genome_ids[organism_id].keys())):
              for contig_id in list(set(genome_ids[organism_id][family])):
                if contig_id not in not_created:
                    not_created.append(contig_id)
                if contig_id in created:
                    created.remove(contig_id)
          else:  
            log.write("'{}' Subtracted file: {}\n".format(organism_id,subtracted_file))
            log.write("'{}' Ready to copy subtracted file to these family folders: {}\n".format(organism_id,", ".join(sorted(list(set(genome_ids[organism_id].keys()))))))
            for family in sorted(list(set(genome_ids[organism_id].keys()))):
                if not family=="":
                    fasta_folder=os.path.join(args.o,family)
                else:
                    fasta_folder=args.o
                if not os.path.exists(os.path.join(fasta_folder, "without_element")):
                    os.makedirs(os.path.join(fasta_folder, "without_element"))
                new_file=os.path.join(fasta_folder, "without_element", "{}_without_element.fasta".format(organism_id))
                unique_contig_ids = set(contig_id for organism_id in genome_ids for contig_id in genome_ids[organism_id].values() for family in genome_ids[organism_id].values() for contig_id in family)
                if os.path.isfile(subtracted_file) and not os.path.isfile(new_file):
                    copy_fasta(subtracted_file,new_file, "subtracted",",".join(unique_contig_ids))
            for contig_id in list(set(sseq.keys())):
              if contig_id not in created:
                  created.append(contig_id)
              if contig_id in not_created:
                  not_created.remove(contig_id)
family="" 
for family in list(set(contig_list.keys())):
  log.close()
  log=open(os.path.join(args.o,"file.log"), 'a')
  log.write("\nFAMILY: {}\n".format(family.upper())) 
  if not family=="":
    fasta_folder=os.path.join(args.o,family)
  else:
    fasta_folder=args.o
  for contig_id in list(set(contig_list[family])):
    besthit=dict()    
    isconta_base=False
    log.write("\nCONTIG_ID: {}\n".format(contig_id.upper()))
    if not contig_id in hits_by_qseqid.keys():
      log.write("'{}' ATTENTION: Element did not have any hits valid (pident>={}% and length>={} bp) for extraction!\n".format(contig_id,str(98),str(3000)))
    else:
      for organism_id in hits_by_qseqid[contig_id].keys():
        log.write("ORGANISM_ID: {}\n".format(organism_id.upper()))
        contig_size=get_value(old, "_".join(contig_id.split("_")[:-1]),"contig size","contig id")
        contig_id_by_size[contig_id]=conta_bases(genomes_files[organism_id], contig_size,contig_id,"contig_id")
        if contig_id_by_size[contig_id]==None:
          if not "ATTENTION: New contig_id not found based on size" in ERROR.keys():
            ERROR["ATTENTION: New contig_id not found based on size"]=[]
          ERROR["ATTENTION: New contig_id not found based on size"].append(contig_id)
          log.write("'{}' ATTENTION: New contig_id not found based on size!\n".format(contig_id)) 
        else:
          if "chromosome" in  contig_id_by_size[contig_id]:
            contig_id_by_size[contig_id]=contig_id_by_size[contig_id].replace("chromosome","")
          with open(os.path.join(args.o, "selected_new_names.tab"), 'r+') as f:
            lines = f.readlines()
            found = False
            for i, line in enumerate(lines):
                columns = line.replace("\n", "").split("\t")
                if "Family" in line and "Old" in line and "New" in line and "Organism" in line and "Method" in line:
                    lines.remove(line)
                if family in line and contig_id in line and contig_id_by_size[contig_id] in line and organism_id in line:
                    if not "Size" in line:
                        columns[0] = family
                        columns[1] = contig_id
                        columns[2] = contig_id_by_size[contig_id]
                        columns[3] = organism_id
                        columns[4] = organism_names[organism_id]
                        columns[5] = columns[5].replace(" ","").split(",")
                        columns[5].append("Size")
                        columns[5] = ",".join(sorted(list(set(columns[5]))))
                        lines[i] = "\t".join(columns) + '\n'
                    found = True
                    break
            if not found:
                lines.append("{}\t{}\t{}\t{}\t{}\tSize\n".format(family, contig_id, contig_id_by_size[contig_id], organism_id,organism_names[organism_id]))
            max_lengths = [max(len(column) for column in [line.split("\t")[i] for line in lines]) for i in range(6)]
            f.seek(0)
            f.write("{:<{}}\t{:<{}}\t{:<{}}\t{:<{}}\t{:<{}}\t{:<{}}\n".format("Family", max_lengths[0], "Old", max_lengths[1], "New", max_lengths[2], "Organism ID", max_lengths[3], "Organism name", max_lengths[4], "Method", max_lengths[5]))
            for line in sorted(list(set(lines))):
                columns = line.replace("\n", "").split("\t")
                f.write("{:<{}}\t{:<{}}\t{:<{}}\t{:<{}}\t{:<{}}\t{:<{}}\n".format(columns[0], max_lengths[0], columns[1], max_lengths[1], columns[2], max_lengths[2], columns[3], max_lengths[3], columns[4], max_lengths[4], columns[5], max_lengths[5]))
            f.truncate()
          if not contig_id==contig_id_by_size[contig_id]:
            log.write("'{}' New contig_id based on size: {}\n".format(contig_id,contig_id_by_size[contig_id]))
        for hit in hits_by_qseqid[contig_id][organism_id]:
          if float(float(hit["pident"]))>=98 and int(int(hit["length"]))>=3000:
            log.write("'{}' BLASTN hit:{}\n".format(contig_id, hit))
            element_start=dict()
            element_end=dict()
            names_list=dict()
            sstrand=dict()
            sstart = int(hit["sstart"])
            send = int(hit["send"])
            if contig_id not in element_start.keys():
                element_start[contig_id] = []
            element_start[contig_id] = sstart
            if contig_id not in element_end.keys():
                element_end[contig_id] = []
            element_end[contig_id] = send
            sseqid = hit["sseqid"]
            log.write("'{}' New contig_id based on BLASTN: {}\n".format(contig_id, sseqid))
            if hit is not None and sseqid != contig_id_by_size[contig_id] and contig_id_by_size[contig_id] is not None:
                log.write("'{}' ATTENTION: Contig_id found by size and by BLASTN search was not a match!\n".format(contig_id))
                if "ATTENTION: Contig_id found by size and by BLASTN search was not a match" not in ERROR.keys():
                    ERROR["ATTENTION: Contig_id found by size and by BLASTN search was not a match"] = []
                ERROR["ATTENTION: Contig_id found by size and by BLASTN search was not a match"].append(contig_id)
            slen = int(hit["slen"])
            qlen = int(hit["qlen"])
            pident = float(hit["pident"])
            align_length = int(hit["length"])
            sstrand[contig_id] = str(hit["sstrand"])
            sstrand[contig_id].lower()
            perc_align_length = int(hit["length"])/int(hit["qlen"])*100
            if qlen is not None:
                if int(align_length) != int(qlen):
                    log.write("'{}' ATTENTION: Alignment length ({}) is not equal to element size({})!\n".format(contig_id, str(align_length), str(qlen)))
            if "chromosome" in sseqid:
                sseqid = sseqid.replace("chromosome", "")
            cid_list.append(sseqid)
            names_list[contig_id] = sseqid
            log.write("'{}' Coordinates of the element based on BLASTN: {}-{}\n".format(contig_id, str(sstart), str(send)))
            family=next((fam for fam, contigs in contig_list.items() if contig_id in contigs), None)
            if not family==None:
              genome_folder=os.path.join(args.o,family,"organism_genomes")
            else:
              genome_folder=os.path.join(args.o,"organism_genomes")
            if contig_id_by_size[contig_id] is None or contig_id_by_size[contig_id] != sseqid:
                new_size = conta_bases(os.path.join(genome_folder, "{}.fasta".format(organism_id)), sseqid, contig_id, "size")
                if new_size is not None:
                    log.write("'{}' New contig size based on BLASTN: {}\n".format(contig_id, new_size))
            fasta_file = os.path.join(genome_folder, organism_id + '.fasta')
            if organism_id in genomes_files.keys():
              exist_fasta=genomes_files[organism_id]
            else:
              exist_fasta=[]
            if organism_id not in genome_ids.keys():
              genome_ids[organism_id]=dict()
            if family not in genome_ids[organism_id]:
              genome_ids[organism_id][family] = []
            if contig_id not in genome_ids[organism_id][family]:
              genome_ids[organism_id][family].append(contig_id)
            if not exist_fasta==[]:
              if os.path.isfile(exist_fasta) and not os.path.isfile(fasta_file):
                copy_fasta(exist_fasta,fasta_file,"organism genome",contig_id)
            if not find_contig_in_file(sseqid, fasta_file, contig_id):
                if "ERROR: Fasta file of sequenced organism genomes does not contain contig_id" not in ERROR.keys():
                    ERROR["ERROR: Fasta file of sequenced organism genomes does not contain contig_id"] = []
                ERROR["ERROR: Fasta file of sequenced organism genomes does not contain contig_id"].append(contig_id)
                if contig_id not in not_created:
                    not_created.append(contig_id)
                log.write("'{}' ERROR: Fasta file of sequenced organism genomes does not contain contig_id!\n".format(contig_id))
                print("'{}' ERROR: Fasta file of sequenced organism genomes does not contain contig_id!".format(contig_id))
            elif (float(align_length)/float(qlen))<0.8:
                log.write("'{}' NOT Contig extraction: Hit alignment length ({}<{}) is lower than 80% of query length ({})!\n".format(contig_id,str(align_length),str(float(0.8*float(qlen))),str(qlen)))
                isbesthit=False
            elif float(pident)<98 or float(align_length)<3000:
                log.write("'{}' NOT Contig extraction: Hit is not valid by pident ({}<{}%) or by alignment length ({}<{})!\n".format(contig_id,str(pident),str(98),str(align_length),str(3000)))
                isbesthit=False
            elif float(pident)>=98 and float(align_length)>=3000 and float(perc_align_length)>=98:
                if not os.path.isdir(os.path.join(fasta_folder, "element")):
                    os.makedirs(os.path.join(fasta_folder, "element"))
                isbesthit=False
                if contig_id not in besthit.keys() or not os.path.isfile(os.path.join(fasta_folder,"element", contig_id + '.fasta')):
                  besthit[contig_id]={"pident":float(pident),"length":int(align_length)}
                  isbesthit=True
                  log.write("'{}' Contig extraction: No hits were extracted for this contig, this hit will be extracted!\n".format(contig_id))
                else:
                  old_pident=float(besthit[contig_id]["pident"])
                  old_align_length=int(besthit[contig_id]["length"])
                  if pident<=old_pident and align_length<=old_align_length:
                    log.write("'{}' NOT Contig extraction: A hit with a larger pident ({}<{}) or larger alignment size ({}<{}) has already been extracted, this hit will NOT be extracted!\n".format(contig_id,str(pident),str(old_pident),str(align_length),str(old_align_length)))
                    isbesthit=False
                  elif pident>=old_pident and align_length>=old_align_length:
                    isbesthit=True
                    besthit[contig_id]={"pident":float(pident),"length":int(align_length)}
                    log.write("'{}' Contig extraction: This hit has a larger pident ({}>{}) or larger alignment length ({}>{}) than the extracted hit, this hit will be extracted!\n".format(contig_id,str(pident),str(old_pident),str(align_length),str(old_align_length)))
                if isbesthit:
                  if not extract_element(organism_id, fasta_folder, sseqid, contig_id, sstart, send, sstrand[contig_id]):
                      if "ERROR: Exception when extracting contig" not in ERROR.keys():
                          ERROR["ERROR: Exception when extracting contig"] = []
                      ERROR["ERROR: Exception when extracting contig"].append(contig_id)
                      log.write("'{}' ERROR: Exception when extracting contig!\n".format(contig_id))
                      print("'{}' ERROR: Exception when extracting contig!".format(contig_id))
                      if contig_id in created:
                          created.remove(contig_id)
                      if contig_id not in not_created:
                          not_created.append(contig_id)
                  elif not os.path.exists(os.path.join(fasta_folder, "element", contig_id + '.fasta')):
                      log.write("'{}' ERROR: Fasta file of the extracted contig was not created!\n".format(contig_id))
                      print("'{}' ERROR: Fasta file of the extracted contig was not created!".format(contig_id))
                      if "ERROR: Fasta file of the extracted contig was not created" not in ERROR.keys():
                          ERROR["ERROR: Fasta file of the extracted contig was not created"] = []
                      ERROR["ERROR: Fasta file of the extracted contig was not created"].append(contig_id)
                      if contig_id in created:
                          created.remove(contig_id)
                      if contig_id not in not_created:
                          not_created.append(contig_id)
                  elif not os.path.getsize(os.path.join(fasta_folder, "element", contig_id + '.fasta')) > 0:
                      log.write("'{}' ERROR: Fasta file of the extracted contig is empty!\n".format(contig_id))
                      print("'{}' ERROR: Fasta file of the extracted contig is empty!".format(contig_id))
                      if "ERROR: Fasta file of the extracted contig is empty" not in ERROR.keys():
                          ERROR["ERROR: Fasta file of the extracted contig is empty"] = []
                      ERROR["ERROR: Fasta file of the extracted contig is empty"].append(contig_id)
                      if contig_id in created:
                          created.remove(contig_id)
                      if contig_id not in not_created:
                          not_created.append(contig_id)
                  else:
                      log.write("'{}' Extracted file: {}\n".format(contig_id, os.path.join(fasta_folder,"element", contig_id + '.fasta')))
                      if contig_id not in created:
                          created.append(contig_id)
                      if contig_id in not_created:
                          not_created.remove(contig_id)    
log=open(os.path.join(args.o,"file.log"), 'a')
created=list(set(created))
if not len(created)==0:
    log.write("\n{} FASTA files were created successfully:{}\n".format(len(created),sorted(created)))
    print("{} FASTA files were created successfully!".format(len(created)))
for key in list(set(ERROR.keys())):
  ERROR[key]=list(set(ERROR[key]))
  log.write("\n{}: {} {}: {}\n".format(key.split(": ")[0],len(ERROR[key]),key.split(": ")[-1],sorted(ERROR[key])))
  print("{}: {} {}!".format(key.split(": ")[0],len(ERROR[key]),key.split(": ")[-1]))
#not_created=list(set(not_created))
if not not_created==[]:
  if len(not_created)==1:
    print("\nERROR: {} Fasta file was not created successfully!".format(len(not_created),not_created))
    log.write("\nERROR: {} Fasta file was not created successfully: {}\n".format(len(not_created),sorted(not_created)))
  else:
    print("\nERROR: {} Fasta files were not created successfully!".format(len(not_created)))
    log.write("\nERROR: {} Fasta files were not created successfully: {}\n".format(len(not_created),sorted(not_created)))
contigs = []
for sublist in contig_list.values():
    contigs.extend(sublist)
contigs=list(set(contigs))
if not len(contigs)==0:
    log.write("\n{} contigs were processed: {}\n".format(len(contigs),sorted(contigs)))
    print("{} contigs were processed!".format(len(contigs)))
for family in sorted(list(set(genomes.keys()))):
  genomes_list=list(set(genomes[family].keys()))
  if "Not found" in genomes_list:
    ngenomes=len(genomes_list)-1
  else:
    ngenomes=len(genomes_list)
  log.write("\n{} organism genomes and {} elements of '{}':{}\n".format(str(ngenomes),str(len(contig_list[family])),family,genomes[family]))
execution=datetime.datetime.now() - start_time
print("\nExecution time: {}".format(execution))
log.write("\nExecution time: {}\n".format(execution))
print("\nExecution time per file: {}".format(execution/len(contigs)))
log.write("\nExecution time per file: {}\n".format(execution/len(contigs)))
new_names.close()
#print("DATETIME: gres.py: {}\n".format(datetime.datetime.now()))
table=[]
header = ["Family", "Input fasta files", "FASTA elementos", "organisms", "without element", "Comentários"]
table.append(header)
for family in contig_list.keys():
    row = [family]
    #element_files_dir = os.path.join(args.e, "Family_" + family, "element/")
    #element_files = [f.name.replace("_element.fasta", "") for f in os.scandir(element_files_dir) if f.name.endswith("_element.fasta")]
    element_files=contig_list[family]
    row.append(len(element_files))
    fasta_elements_dir = os.path.join(args.o, family, "element/")
    fasta_elements = [f.name.replace(".fasta", "") for f in os.scandir(fasta_elements_dir) if f.name.endswith(".fasta")]
    row.append(len(fasta_elements))
    organisms_dir = os.path.join(args.o, family, "organism_genomes/")
    organisms = [f.name.replace(".fasta", "") for f in os.scandir(organisms_dir) if f.name.endswith(".fasta")]
    row.append(len(organisms))
    without_element_dir = os.path.join(args.o, family, "without_element/")
    if not os.path.isdir(without_element_dir):
       os.makedirs(without_element_dir)
    without_element = [f.name.replace("_without_element.fasta", "") for f in os.scandir(without_element_dir) if f.name.endswith("_without_element.fasta")]
    row.append(len(without_element))
    comments = []
    if len(element_files)!= len(fasta_elements):
        only_in_element_files = set(element_files) - set(fasta_elements)
        only_in_fasta_elements = set(fasta_elements) - set(element_files)
        if only_in_element_files==set() and not only_in_fasta_elements==set():
          comments.append("only in fasta elements: {}".format(only_in_fasta_elements))
        elif only_in_fasta_elements==set() and not element_files==set():
          comments.append("only in input fasta files:{}".format(only_in_element_files))
        else:
          comments.append("only in input fasta files:{}, only in fasta elements:{}".format(only_in_element_files,only_in_fasta_elements))
    if len(organisms)!= len(without_element):
        only_in_organisms = set(organisms) - set(without_element)
        only_in_without_element = set(without_element) - set(organisms)
        if only_in_organisms==set() and not only_in_without_element==set():
          comments.append("only in without element:{}".format(only_in_without_element))
        elif only_in_without_element==set() and not organisms==set():
          comments.append("only in organisms:{}".format(only_in_organisms))
        else:
          comments.append("only in organisms:{}, only in without element:{}".format(only_in_organisms,only_in_without_element))
    row.append(", ".join(comments))
    table.append(row)
max_lengths = [max(len(str(row[i])) for row in table) for i in range(len(table[0]))]
for row in table:
    for i, cell in enumerate(row):
        print(str(cell).ljust(max_lengths[i] + 2), end="")
        log.write(str(cell).ljust(max_lengths[i] + 2))
    print()
    log.write("\n")
print("\nEnd")
log.write("\nEnd\n")
log.close()
exit()