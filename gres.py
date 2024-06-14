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

version="1.3.2"
print('gres v{}'.format(version))

start_time = datetime.datetime.now()
call=os.path.abspath(os.getcwd())

contig_list=dict()
element_files=dict()
error=dict()
old=''
new=''
delimiter=''
contig_size=dict()
not_created=[]
created=[]
genome_ids=dict()
cid_list=[]
new_element_id=dict()


parser = argparse.ArgumentParser(
  prog='gres',
  description=
  'Genome Retriever and Element Subtractor',
  epilog='https://github.com/GiulianaPola/gres',add_help=False)
parser.add_argument('-i',help="contig_idation file or folder")
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

def printing_help():
    print('(c) 2023. Giuliana Pola & Arthur Gruber\n')
    print('For more information access: ')
    print('Usage: updateID.py -i <path|file> -old <old>')
    print('Genome Retriever and Element Subtractor')
    print('\nOptional arguments:')
    print('  -h, --help             Show this help message and exit')
    print('  -i <path|file>        contig_idation files or folder (required)')
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
                        delimiter=find_delimiter(lines[0])
                        if delimiter==None:
                          contig_list[""]=lines
                        else:
                          columns=lines[0].replace("\n","").split(delimiter)
                          family_index=""
                          contig_index=""
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
                              if not "" in contig_list.keys():
                                contig_list[""]=[]
                              contig_list[""].append(columns[contig_index])
                              
            else:
                print("Contig_id's list (-i) is empty.")
                exit()
        else:
            print("Contig_id's list (-i) does not exist.")
            exit()
    except Exception as e:
        print("Error Contig_id's list (-i): ", str(e))
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
           print("{}: {}\n".format(i, line))
        exit()
    
    args.e=get_actual_path(args.e,call)
    if os.path.isdir(args.e):
      if not args.e[-1]=="/":
        args.e=args.e+"/"
      if not os.path.exists(args.e):
          print("ERROR: Element's fasta files' folder '{}' does not exist!".format(args.e))
          exit()
      elif not os.listdir(args.e):
        print("ERROR: Element's fasta files' folder '{}' is empty!".format(args.e))
        exit()
      elif os.path.isdir(args.e):
        for root, dirs, files in os.walk(args.e):
          for file in files:
              if file.endswith('.fasta'):
                  file=os.path.join(root, file)
                  if os.path.split(file)[1] not in element_files:
                    element_files[os.path.split(file)[1]]=file
        if element_files==[]:
          print("ERROR: Input folder '{}' does not contain element's fasta files (.fasta)!".format(args.e))
          exit()
    elif os.path.isfile(args.e):
      if not os.path.exists(args.e):
          print("ERROR: Element's fasta file '{}' does not exist!".format(args.e))
          exit()
      element_files[os.path.split(args.e)[1]]=args.e
    
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
      print(e)
      print("ERROR: Output folder was not written!")
      exit()
    return args.i,old,new
               
def get_value(table, row,wanted_col, search_col):
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
            if "ERROR: get_value exception" not in error.keys():
              error['ERROR get_value exception']=[]
            error['ERROR get_value exception'].append(contig_id)              
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

def download_sequences(organism_id, contig_id):
    cmd = "p3-genome-fasta {} > {}.fasta".format(organism_id,organism_id)
    log.write("'{}' Download the fasta file of sequenced organism genomes: {}\n".format(contig_id,cmd))
    
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        log.write("ERROR '{}': Failed to download organism_id '{}' fasta due to {}!\n".format(contig_id,organism_id, str(e.output)))
        return False
    except Exception as e:
        log.write("ERROR '{}': Failed to download organism_id '{}' fasta due to {}!\n".format(contig_id,organism_id, e))
        return False
    else:
        log.write("'{}' Fasta from organism_id '{}' was downloaded!\n".format(contig_id,organism_id))
        return True
        


def conta_bases(fasta_file, filtering,contig_id,wanted):
    try:
        cmd = "conta_bases.pl -i '{}' |grep '{}'".format(fasta_file, filtering)
        
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        output, error = process.communicate()
    
        if error or not output:
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
        if "ERROR: conta_bases.pl exception" not in error.keys():
          error['ERROR: conta_bases.pl exception']=[]
        error['ERROR: conta_bases.pl exception'].append(contig_id)              
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
            print("{}: {}".format(i, line))

def find_contig_in_file(new_contig_id, ffile,contig_id):
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
    try:
        count = int(subprocess.check_output("grep -c '^>' " + filename, shell=True))
        return count
    except Exception as e:
        print("'{}' ERROR: conta_bases.pl exception: {}".format(contig_id,e))
        if "ERROR: conta_bases.pl exception" not in error.keys():
          error['ERROR: conta_bases.pl exception']=[]
        error['ERROR: conta_bases.pl exception'].append(contig_id)              
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
            print("{}: {}".format(i, line))
        return None

def get_sequences(family,contig_id, fasta_folder, old,new):
    new_names=open(os.path.join(args.o,"selected_new_names.tab"), 'r+')
    coordinates = []
    sequence = ''
    organism_name=None
    try:
        fasta_folder=os.path.abspath(fasta_folder)
        log.write("\nCONTIG_ID: {}\n".format(contig_id))
        log.write("'{}' fasta folder: {}\n".format(contig_id,fasta_folder))
        if os.path.split(args.old)[-1]=="":
          tail="table -old"
        else:
          tail=os.path.split(args.old)[-1]+" (-old)"
        if not find_contig_in_file(contig_id.split("_")[0], args.old,contig_id):
          log.write("'{}' ERROR: Contig_id not found in table (-old)!\n".format(contig_id,tail))
          if not "ERROR: Contig_id not found in {}".format(tail) in error.keys():
              error["ERROR: Contig_id not found in {}".format(tail)]=[]
          error["ERROR: Contig_id not found in {}".format(tail)].append(contig_id)
        else:
          organism_name=get_value(old, contig_id.split("_")[0],"Organism name", "contig id")
          if organism_name==None:
            log.write("ERROR '{}': organism_name was not found in {}!\n".format(contig_id,tail))
            if not "ERROR: organism_name was not found in {}".format(tail) in error.keys():
              error["ERROR: organism_name was not found".format(tail)]=[]
            error["ERROR: organism_name was not found".format(tail)].append(contig_id)
          else:
            log.write("'{}' Organism_name ({}): {}\n".format(contig_id,tail,organism_name))
            values_old=dict()
            element_size=get_value(old, contig_id.split("_")[0],"Element size","contig id")
            contig_size=get_value(old, contig_id.split("_")[0],"contig size","contig id")
            old_organism_id=get_value(old, contig_id.split("_")[0],"Organism id", "contig id")
            if not element_size==None:
              values_old["element_size"]=element_size
            if contig_size==None:
              log.write("ERROR '{}': contig_size was not found in {}!\n".format(contig_id,tail))
              if not "ERROR: contig_size was not found in {}".format(tail) in error.keys():
                error["ERROR: contig_size was not found in {}".format(tail)]=[]
              error["ERROR: contig_size was not found in {}".format(tail)].append(contig_id)
            else:
              values_old["contig_size"]=contig_size 
              if not old_organism_id==None:
                log.write("'{}' Organism_id ({}): {}\n".format(contig_id,tail,old_organism_id))
                values_old["organism_id"]=old_organism_id
            if len(values_old.keys())>0:
              log.write("'{}' {}: {}\n".format(contig_id,tail,values_old))
            if not organism_name==None:
              organism_id=get_value(new, organism_name,"genome id","genome name")
              if os.path.split(args.old)[-1]=="":
                tail="table -new"
              else:
                tail=os.path.split(args.new)[-1]+" (-new)"
            if not organism_id==None and not old_organism_id==None:
              if not old_organism_id==organism_id:
                log.write("'{}' ATTENTION: organism_id has changed!\n".format(contig_id))
            if organism_id==None:
              log.write("'{}' find organism_id in file: grep '{}' '{}'\n".format(contig_id, organism_name, args.new))
              log.write("ERROR '{}': Organism_id was not found in {}!\n".format(contig_id, tail))
              if not "ERROR: Organism_id was not found in {}".format(tail) in error.keys():
                error["ERROR: Organism_id was not found in {}".format(tail)]=[]
              error["ERROR: Organism_id was not found in {}".format(tail)].append(contig_id)
            else:
              if os.path.split(args.new)[-1]=="":
                tail="-new"
              else:
                tail=os.path.split(args.new)[-1]+" (-new)"
              log.write("'{}' New organism_id ({}): {}\n".format(contig_id,tail,organism_id))
              if not os.path.exists(fasta_folder.replace("\ "," ")):
                os.makedirs(fasta_folder.replace("\ "," "))
              fasta_folder=fasta_folder.replace("\\ "," ")
              genome_folder=os.path.join(fasta_folder,"organism_genomes")
              if not os.path.exists(genome_folder):
                  os.makedirs(genome_folder)
              os.chdir(genome_folder)
              if not download_sequences(organism_id,contig_id):  
                if not "ERROR: Sequenced genomes of the organism have not been downloaded" in error.keys():
                  error["ERROR: Sequenced genomes of the organism have not been downloaded"]=[]
                error["ERROR: Sequenced genomes of the organism have not been downloaded"].append(contig_id)
                log.write("'{}' ERROR: Sequenced genomes of the organism have not been downloaded!\n".format(contig_id))
              elif not os.path.exists(os.path.join(genome_folder,"{}.fasta".format(organism_id))):
                if not "ERROR: Fasta file of sequenced genomes of the organism was not created" in error.keys():
                  error["ERROR: Fasta file of sequenced genomes of the organism was not created"]=[]
                error["ERROR: Fasta file of sequenced genomes of the organism was not created"].append(contig_id)
                log.write("'{}' ERROR: Fasta file of sequenced genomes of the organism was not created!\n".format(contig_id))
              elif not os.path.getsize(os.path.join(genome_folder,"{}.fasta".format(organism_id))) > 0:
                if not "ERROR: Fasta file of sequenced genomes of the organism is empty" in error.keys():
                  error["ERROR: Fasta file of sequenced genomes of the organism is empty"]=[]
                error["ERROR: Fasta file of sequenced genomes of the organism is empty"].append(contig_id)
                log.write("'{}' ERROR: Fasta file of sequenced genomes of the organism is empty!\n".format(contig_id))
              else:
                log.write("'{}' Organism's genome fasta: {}\n".format(contig_id,os.path.join(genome_folder,"{}.fasta".format(organism_id))))
                cid_list=[]
                if not contig_size==None:
                  new_contig_id=conta_bases(os.path.join(genome_folder,"{}.fasta".format(organism_id)), contig_size,contig_id,"contig_id")                         
                  if new_contig_id==None:
                    if not "ATTENTION: New contig_id not found based on size" in error.keys():
                      error["ATTENTION: New contig_id not found based on size"]=[]
                    error["ATTENTION: New contig_id not found based on size"].append(contig_id)
                    log.write("'{}' ATTENTION: New contig_id not found based on size!\n".format(contig_id)) 
                  else:
                    if "chromosome" in  new_contig_id:
                      new_contig_id=new_contig_id.replace("chromosome","")
                    if not "{}\t{}\t{}\n".format(family,contig_id,new_contig_id) in new_names.read():
                      new_names.write("{}\t{}\t{}\n".format(family,contig_id,new_contig_id))
                      new_names.close()
                      new_names=open(os.path.join(args.o,"selected_new_names.tab"), 'r+')
                    if not contig_id==new_contig_id:
                      log.write("'{}' New contig_id based on size: {}\n".format(contig_id,new_contig_id)) 
                    cid_list.append(new_contig_id)
                filename=''
                namesize=9999
                positive=None
                for name in element_files.keys():
                  if contig_id.lower() in name.lower():
                    if len(name)<namesize:
                      filename=name
                      namesize=len(name)
                if filename=='':
                  if os.path.basename(os.path.normpath(args.e))=="":
                    tail="fasta folder (-e)"
                  else:
                    tail=os.path.basename(os.path.normpath(args.e))+" (-e)"
                    
                  log.write("'{}' Find element's fasta in folder: find {} -type f -name '*{}*.fasta'\n".format(contig_id,args.e,contig_id))
                  log.write("'{}' ERROR: Element's fasta was not found in {}!\n".format(contig_id,tail))
                  if not "ERROR: Element's fasta was not found in {}".format(tail) in error.keys():
                    error["ERROR: Element's fasta was not found in {}".format(tail)]=[]
                  error["ERROR: Element's fasta was not found in {}".format(tail)].append(contig_id)
                else:
                  if not os.path.exists(os.path.join("/".join(fasta_folder.split("/")[0:-1]),'BLASTN').replace("\ "," ")):
                    os.mkdir(os.path.join("/".join(fasta_folder.split("/")[0:-1]),'BLASTN').replace("\ "," "))
                  family_folder = os.path.join("/".join(fasta_folder.split("/")[0:-1])).replace("\ "," ")
                  positive=blastn_element(family_folder, element_files[filename], os.path.join(genome_folder,"{}.fasta".format(organism_id)),contig_id)
                  if positive==None:
                    if not "ERROR: No hits found" in error.keys():
                      error["ERROR: No hits found"]=[]
                    error["ERROR: No hits found"].append(contig_id)
                  else:
                    sstart=positive["organism genome"]["sstart"]
                    send=positive["organism genome"]["send"]
                    if not contig_id in element_start.keys():
                        element_start[contig_id]=[]
                    element_start[contig_id]=sstart
                    if not contig_id in element_end.keys():
                        element_end[contig_id]=[]
                    element_end[contig_id]=send
                    sseqid=positive["organism genome"]["sseqid"]
                    if not positive==None and not sseqid==new_contig_id and not new_contig_id==None:
                      log.write("'{}' ERROR: Contig_id found by size and by BLASTN search was not a match!\n".format(contig_id))
                      if not "ERROR: Contig_id found by size and by BLASTN search was not a match" in error.keys():
                        error["ERROR: Contig_id found by size and by BLASTN search was not a match"]=[]
                      error["ERROR: Contig_id found by size and by BLASTN search was not a match"].append(contig_id)
                    else:
                      align_length=positive["alignment"]["length"]
                      element_seq=positive["organism genome"].pop("sseq")
                      #sseq[sseqid].append(element_seq)
                      log.write("'{}' Blast result:{}\n".format(contig_id,positive))
                      slen=positive["organism genome"]["slen"]
                      qlen=positive["element"]["qlen"]
                      sstrand[contig_id]=str(positive["organism genome"]["sstrand"])
                      sstrand[contig_id].lower()
                      if not element_size==None:
                        if not int(align_length)==int(qlen):
                          log.write("'{}' ATTENTION: Alignment length ({}) is not equal to element size({})!\n".format(contig_id, str(align_length),str(qlen)))
                      if "chromosome" in positive:
                        sseqid=sseqid.replace("chromosome","")
                      cid_list.append(sseqid)
                      if not "{}\t{}\t{}\n".format(family,contig_id,sseqid) in new_names.read():
                        new_names.write("{}\t{}\t{}\n".format(family,contig_id,sseqid))
                      names_list[contig_id]=sseqid
                      log.write("'{}' New contig_id based on BLASTN: {}\n".format(contig_id,sseqid))
                      log.write("'{}' Coordinates of the element based on BLASTN: {}-{}\n".format(contig_id,str(sstart),str(send)))
                      if new_contig_id==None or not new_contig_id==sseqid:
                          new_size=conta_bases(os.path.join(genome_folder,"{}.fasta".format(organism_id)), sseqid,contig_id,"size")
                          if not new_size==None:
                            log.write("'{}' New contig size based on BLASTN: {}\n".format(contig_id,new_size))
                    
                      for new_cid in list(dict.fromkeys(cid_list)):
                        if all(item is None for item in list(dict.fromkeys(cid_list))):
                          if contig_id not in not_created:
                            not_created.append(contig_id)
                          if contig_id in created:
                            created.remove(contig_id)
                        if not new_cid==None:
                          if "chromosome" in new_cid:
                            new_cid=new_cid.replace("chromosome","")
                          if not organism_id in new_element_id.keys():
                            new_element_id[organism_id]=[]
                          new_element_id[organism_id].append(new_cid)
                          fasta_file = os.path.join(genome_folder, organism_id + '.fasta')
                          if not find_contig_in_file(new_cid, fasta_file,contig_id):
                            if not "ERROR: Fasta file of sequenced organism genomes does not contain contig_id" in error.keys():
                              error["ERROR: Fasta file of sequenced organism genomes does not contain contig_id"]=[]
                            error["ERROR: Fasta file of sequenced organism genomes does not contain contig_id"].append(contig_id)
                            if contig_id not in not_created:
                              not_created.append(contig_id)
                            log.write("'{}' ERROR: Fasta file of sequenced organism genomes does not contain contig_id!\n".format(contig_id))
                          else:
                            if not os.path.isdir(os.path.join(fasta_folder, "element")):
                              os.makedirs(os.path.join(fasta_folder, "element"))
                            if not extract_element(organism_id, fasta_folder, new_cid,contig_id,sstart,send,sstrand[contig_id]):
                              if not "ERROR: Exception when extracting contig" in error.keys():
                                error["ERROR: Exception when extracting contig"]=[]
                              error["ERROR: Exception when extracting contig"].append(contig_id)
                              log.write("'{}' ERROR: Exception when extracting contig!\n".format(contig_id))  
                              if contig_id in created:
                                  created.remove(contig_id)
                              if contig_id not in not_created:
                                not_created.append(contig_id)
                            elif not os.path.exists(os.path.join(fasta_folder,"element", contig_id + '.fasta')):
                              log.write("'{}' Extracted file: {}\n".format(contig_id,os.path.join(fasta_folder, contig_id + '.fasta')))
                              log.write("'{}' ERROR: Fasta file of the extracted contig was not created!\n".format(contig_id))  
                              if not "ERROR: Fasta file of the extracted contig was not created" in error.keys():
                                error["ERROR: Fasta file of the extracted contig was not created"]=[]
                              error["ERROR: Fasta file of the extracted contig was not created"].append(contig_id)
                              if contig_id in created:
                                  created.remove(contig_id)
                              if contig_id not in not_created:
                                not_created.append(contig_id)
                            elif not os.path.getsize(os.path.join(fasta_folder,"element", contig_id + '.fasta')) > 0:
                              log.write("'{}' ERROR: Fasta file of the extracted contig is empty!\n".format(contig_id))  
                              if not "ERROR: Fasta file of the extracted contig is empty" in error.keys():
                                error["ERROR: Fasta file of the extracted contig is empty"]=[]
                              error["ERROR: Fasta file of the extracted contig is empty"].append(contig_id)
                              if contig_id in created:
                                  created.remove(contig_id)
                              if contig_id not in not_created:
                                not_created.append(contig_id)
                            else:
                              if contig_id not in created:
                                  created.append(contig_id)
                              if contig_id in not_created:
                                not_created.remove(contig_id)
                          
    except Exception as e:
        if "gres Exception" not in error.keys():
          error["gres Exception"]=[]
        error["gres Exception"].append(contig_id)
        print("ERROR '{}': {}".format(contig_id, e))
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
            print("{}: {}".format(i, line))
    os.chdir(call)
    new_names.close()

def run_blastn(contig_id, element_fasta, organism_multifasta, outfile, outfmt_str):
    try:
        blastn_cline = NcbiblastnCommandline(query=element_fasta, subject=organism_multifasta, evalue=0.05, num_alignments=100, reward=2, penalty=-3, gapopen=5, gapextend=2, dust="yes", out=outfile, outfmt='"6 {}"'.format(outfmt_str))
        log.write("'{}' Running BLASTN command: {}\n".format(contig_id,blastn_cline))
        stdout, stderr = blastn_cline()
    except Exception as e:
        print("'{}' ERROR BLASTN: {}".format(contig_id, e))
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
            print("{}: {}".format(i, line))
        return False
    else:
        return True

def blastn_element(family_folder, element_fasta, organism_multifasta, contig_id):
    outfile = os.path.join(family_folder, 'BLASTN', 'blastn_{}.txt'.format(contig_id))
    result = dict()
    outfmt_str = "qseqid qlen qstart qend sseqid sstrand slen sstart send pident length mismatch gapopen evalue bitscore sseq"
    cols = outfmt_str.split()
    try:
        os.chdir(os.path.join(family_folder, 'BLASTN'))
        if run_blastn(contig_id, element_fasta, organism_multifasta, outfile, outfmt_str) and os.path.isfile(outfile):
            if os.path.getsize(outfile) > 0:
                with open(outfile, 'r') as f:
                    hit = f.readline()
                result = {'element': {}, 'organism genome': {}, 'alignment': {}}
                for col, value in zip(cols, hit.split()):
                    if col.startswith('q'):
                        result['element'][col] = value
                    elif col.startswith('s'):
                        result['organism genome'][col] = value
                    else:
                        result['alignment'][col] = value
                return result
            else:
                log.write("'{}' ERROR BLASTN: No hits found!\n".format(contig_id))
                return None
        else:
            log.write("'{}' ERROR BLASTN: Result file was not created!\n".format(contig_id))
            return None
    except Exception as e:
        print("'{}' ERROR BLASTN: {}".format(contig_id, e))
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
            print("{}: {}".format(i, line))
        return None

def extract_element(organism_id, fasta_folder, new_contig_id,old_contig_id,start,end,sstrand):
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
        log.write("'{}' Extracted file: {}\n".format(old_contig_id,output_file))
        return True
      else:
        log.write("'{}' ERROR: Contig {} not found in FASTA file!\n".format(old_contig_id,new_contig_id))
        return False
    except Exception as e: 
        print("'{}' ERROR: {}\n".format(old_contig_id, str(e)))
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
          print("{}: {}".format(i, line))
        return False
    else:
        return True

                    
      
def write_records_to_fasta(records, output_file):
    try:
        SeqIO.write(records, output_file, "fasta")
    except Exception as e: 
      print("'{}' ERROR: {}\n".format(organism_id, str(e)))
      exc_type, exc_value, exc_traceback = sys.exc_info()
      formatted_lines = traceback.format_exc().splitlines()
      for i, line in enumerate(formatted_lines):
          print("{}: {}".format(i, line))
    else:
      log.write("'{}' Subtracted file: {}\n".format(organism_id,output_file))

def subtract_element(organism_id, fasta_folder, names_list, start, end, sseq,sstrand):
    try:
        fasta_file = os.path.join(genome_folder, organism_id + '.fasta')
        records = list(SeqIO.parse(fasta_file, "fasta"))
        changed=False

        if not os.path.isdir(os.path.join(fasta_folder, "without_element")):
            os.makedirs(os.path.join(fasta_folder, "without_element"))

        for rec in records:
            rec_old_contig_id = []
            rec_start = [1]
            rec_end = [len(rec.seq)]
            rec_sseq = []
            if rec.id in names_list.values():
                for old_contig_id in names_list.keys():
                    if names_list[old_contig_id] == rec.id:
                        if old_contig_id in sseq.keys():
                            rec_old_contig_id.append(old_contig_id)
                            if int(start[old_contig_id]) <= int(end[old_contig_id]) or sstrand[old_contig_id]=='plus':
                                if int(start[old_contig_id])==1:
                                  rec_end.append(1)
                                else:
                                  rec_end.append(int(start[old_contig_id])-1)
                                rec_start.append(int(end[old_contig_id]))
                            if int(start[old_contig_id]) > int(end[old_contig_id]) or sstrand[old_contig_id]=='minus':
                                if int(end[old_contig_id])==1:
                                  rec_end.append(1)
                                else:
                                  rec_end.append(int(end[old_contig_id])-1)
                                rec_end.append(int(end[old_contig_id])-1)
                                rec_start.append(int(start[old_contig_id]))
                            rec_sseq.append(sseq[old_contig_id])

            if len(rec_start) > 1 and len(rec_end) > 1:
                changed=True
                rec_start = sorted([int(x) for x in rec_start])
                rec_end = sorted([int(y) for y in rec_end])
                log_str =""
                for s, e in zip(rec_start, rec_end):
                  log_str = log_str+"{}-{} ".format(str(s), str(e))
                log.write("'{}' '{}' Subtract element: '{}' {}\n".format(organism_id, ",".join(rec_old_contig_id), rec.id,log_str))
                for contig_id in rec_old_contig_id:
                    rec_seq = re.sub(r'[^ACGTUNacgtun]', '', str(rec.seq).upper())
                    seq = re.sub(r'[^ACGTUNacgtun]', '', str(sseq[contig_id]).upper())
                    if str(rec_seq).find(str(seq))> 0:
                        match=(re.search(str(seq), str(rec_seq)))
                        istart=match.start()
                        iend=match.end()
                        rec.seq = Seq(str(rec_seq).replace(str(seq), ""))
                        if str(rec.seq).find(str(seq))> 0:
                          rec.seq = Seq("".join(str(rec_seq).split(str(seq))))
                          if str(rec.seq).find(str(seq))> 0:
                            seq_list=[]
                            seq_list.append(str(rec_seq)[0:istart])
                            seq_list.append(str(rec_seq)[iend:len(rec_seq)])
                            rec.seq=Seq("".join(seq_list))
                            if str(rec.seq).find(str(seq))> 0:
                              new_seq=""
                              for s,e in zip(rec_start, rec_end):
                                new_seq=new_seq+rec_seq[s-1:e]
                              rec.seq=Seq(new_seq)
                            
                    else:
                        new_seq=""
                        for s,e in zip(rec_start, rec_end):
                          #BLASTN: 675-10932
                          #match 674-10931 vs 1-675 10932-33627
                          #0-674 10931-33627
                          #1-674 10932-33627
                          if s==1 and e==1:
                            new_seq=new_seq+""
                          else:
                            new_seq=new_seq+rec_seq[s-1:e]
                        rec.seq=Seq(new_seq)
                if str(rec.seq).find(str(seq))> 0:
                  log.write("'{}' ERROR: Element sequence '{}' was not subtracted from contig '{}' in organism fasta file!\n".format(organism_id,contig_id,rec.id))
                rec.description=log_str
        output_file = os.path.join(fasta_folder, "without_element", "{}-{}.fasta".format(organism_id, ",".join(sseq.keys())))
        if changed:
          write_records_to_fasta(records, output_file)
    except Exception as e: 
        print("'{}' ERROR: {}\n".format(organism_id, str(e)))
        exc_type, exc_value, exc_traceback = sys.exc_info()
        formatted_lines = traceback.format_exc().splitlines()
        for i, line in enumerate(formatted_lines):
          print("{}: {}".format(i, line))
        return False
    else:
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
    log.write('\nWorking directory: {}\n'.format(call))
    log.write('\nCommand line: {}\n'.format(' '.join(sys.argv)))
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
    new_names.write("Family\tOld\tNew\n")
  except Exception as e:
    print("ERROR: Selected_new_names.tab file was not written!")
    exit()
new_names.close()
genomes=dict()

for family in sorted(list(set(contig_list.keys()))):
    genome_ids=dict()
    log=open(os.path.join(args.o,"file.log"), 'a')
    fasta_folder=""
    genome_folder=""
    log.write("\nFAMILY: {}\n".format(family.upper()))
    os.chdir(call)
    if not family=="":
      fasta_folder=os.path.join(args.o,family,"FASTA")
    else:
      fasta_folder=os.path.join(args.o, "FASTA")
    if not os.path.exists(fasta_folder):
      os.makedirs(fasta_folder.replace("\ "," "))
    contig_list[family]=list(set(contig_list[family]))
    for contig_id in sorted(contig_list[family]):
      log=open(os.path.join(args.o,"file.log"), 'a')
      old_organism_id=get_value(old, contig_id,"Organism id", "contig id")
      organism_name=get_value(old, contig_id.split("_")[0],"Organism name", "contig id")
      organism_id=get_value(new, organism_name,"genome id","genome name")
      if organism_id==None and old_organism_id==None:
          if "Not found" not in genome_ids.keys():
            genome_ids["Not found"]=[]
          if contig_id not in genome_ids.values():
            genome_ids["Not found"].append(contig_id)
            genome_ids["Not found"]=sorted(list(set(genome_ids["Not found"])))
      elif organism_id==None:
          if old_organism_id not in genome_ids.keys():
            genome_ids[old_organism_id]=[]
          if contig_id not in genome_ids.values():
            genome_ids[old_organism_id].append(contig_id)
            genome_ids[old_organism_id]=sorted(list(set(genome_ids[old_organism_id])))
      else:
          if organism_id not in genome_ids.keys():
            genome_ids[organism_id]=[]
          if contig_id not in genome_ids.values():
            genome_ids[organism_id].append(contig_id)
            genome_ids[organism_id]=sorted(list(set(genome_ids[organism_id])))
    log=open(os.path.join(args.o,"file.log"), 'a')
    
    if family not in genomes.keys():
      genomes[family]=genome_ids
    
    for organism_id in sorted(list(set(genome_ids.keys()))):
      log.write("\nORGANISM_ID: {}\n".format(organism_id.upper()))
      element_start=dict()
      element_end=dict()
      new_element_id=dict()
      sseq=dict()
      names_list=dict()
      sstrand=dict()
      log=open(os.path.join(args.o,"file.log"), 'a')
      for contig_id in sorted(list(set(genome_ids[organism_id]))):
        get_sequences(family,contig_id,fasta_folder,old,new)
        log.close()
        log=open(os.path.join(args.o,"file.log"), 'a')
      genome_folder=os.path.join(fasta_folder,"organism_genomes")
      fasta_file = os.path.join(genome_folder, organism_id + '.fasta')
      nseqs=count_sequences(fasta_file)
      if organism_id in new_element_id.keys() and not len(sseq.keys())==0:
        log.write("\n")
        os.chdir(genome_folder)
        if nseqs==None:
          if not "ERROR: Couldn't count the contigs in the organism's multifasta" in error.keys():
            error["ERROR: Couldn't count the contigs in the organism's multifasta"]=[]
          error["ERROR: Couldn't count the contigs in the organism's multifasta"].extend(genome_ids[organism_id])
          log.write("'{}' ERROR: Couldn't count the contigs in the Fasta file of the organism's multifasta!\n".format(organism_id))
        elif nseqs==0:
          if not "ERROR: The organism's multifasta has no contig" in error.keys():
            error["ERROR: The organism's multifasta has no contig"]=[]
          error["ERROR: The organism's multifasta has no contig"].extend(genome_ids[organism_id])
          log.write("'{}' ERROR: The organism's multifasta has no contig!\n".format(organism_id))
          for contig_id in genome_ids[organism_id]:
            if contig_id in created:
              created.remove(contig_id)
            if contig_id not in not_created:
              not_created.append(contig_id)
        else:
          if nseqs==1:
            if not "ATTENTION: The organism's multifasta has only one contig" in error.keys():
              error["ATTENTION: The organism's multifasta has only one contig"]=[]
            error["ATTENTION: The organism's multifasta has only one contig"].append(contig_id)
            log.write("'{}' ATTENTION: The organism's multifasta has only one contig!\n".format(organism_id))
          subtracted_file=os.path.join(fasta_folder, "without_element", "{}-{}.fasta".format(organism_id,",".join(sseq.keys())))
          if not subtract_element(organism_id, fasta_folder, names_list,element_start,element_end,sseq,sstrand):
            if not "ERROR: Exception when subtracting contig" in error.keys():
              error["ERROR: Exception when subtracting contig"]=[]
            error["ERROR: Exception when subtracting contig"].append(organism_id)
            log.write("'{}' ERROR: Exception when subtracting contig!\n".format(organism_id))
            for contig_id in genome_ids[organism_id]:
              if contig_id not in not_created:
                  not_created.append(contig_id)
              if contig_id in created:
                  created.remove(contig_id)
          elif not os.path.exists(subtracted_file):
            log.write("'{}' ERROR: Fasta file of the subtracted contig was not created!\n".format(organism_id,",".join(genome_ids[organism_id])))  
            if not "ERROR: Fasta file of the subtracted contig was not created" in error.keys():
              error["ERROR: Fasta file of the subtracted contig was not created"]=[]
            error["ERROR: Fasta file of the subtracted contig was not created"].extend(genome_ids[organism_id])
            for contig_id in genome_ids[organism_id]:
              if contig_id not in not_created:
                  not_created.append(contig_id)
              if contig_id in created:
                  created.remove(contig_id)
          elif not os.path.getsize(subtracted_file) > 0:
            log.write("'{}' ERROR: Fasta file of the subtracted contig is empty!\n".format(organism_id,",".join(genome_ids[organism_id])))  
            if not "ERROR: Fasta file of the subtracted contig is empty" in error.keys():
              error["ERROR: Fasta file of the subtracted contig is empty"]=[]
            error["ERROR: Fasta file of the subtracted contig is empty"].extend(genome_ids[organism_id])
            for contig_id in genome_ids[organism_id]:
              if contig_id not in not_created:
                  not_created.append(contig_id)
              if contig_id in created:
                  created.remove(contig_id)
          else:  
            for contig_id in sseq.keys():
              if contig_id not in created:
                  created.append(contig_id)
              if contig_id in not_created:
                  not_created.remove(contig_id)
      log.close()
      log=open(os.path.join(args.o,"file.log"), 'a')
created=list(set(created))
if not len(created)==0:
    log.write("\n{} FASTA files were created successfully:{}\n".format(len(created),sorted(created)))
    print("{} FASTA files were created successfully!".format(len(created)))
for key in error.keys():
  log.write("\n{}: {} {}: {}\n".format(key.split(": ")[0],len(error[key]),key.split(": ")[-1],sorted(error[key])))
  print("{}: {} {}!".format(key.split(": ")[0],len(error[key]),key.split(": ")[-1]))
not_created=list(set(not_created))
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
log.write("\n")
for name in genomes.keys():
  genomes_list=list(set(genomes[name].keys()))
  if "Not found" in genomes_list:
    ngenomes=len(genomes_list)-1
  else:
    ngenomes=len(genomes_list)
  log.write("{} organism genomes and {} elements of '{}':{}\n".format(str(ngenomes),str(len(contig_list[name])),name,genomes[name]))
execution=datetime.datetime.now() - start_time
print("\nExecution time: {}".format(execution))
log.write("\nExecution time: {}\n".format(execution))
print("\nExecution time per file: {}".format(execution/len(contigs)))
log.write("\nExecution time per file: {}\n".format(execution/len(contigs)))
print("End")
log.write("\nEnd\n")
log.close()
new_names.close()
exit()