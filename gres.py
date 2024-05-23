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

version="1.3.0"
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
contigs=[]
created=[]
genome_ids=dict()

try:
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
        else:
            if organism_id not in organism.keys():
              genome_ids[organism_id]=[]
            genome_ids[organism_id].append(contig_id)
            log.write("'{}': Fasta from organism_id '{}' was downloaded!\n".format(contig_id,organism_id))
            return True
            


    def count_contig_bases(fasta_file, filtering,contig_id,wanted):
        try:
            cmd = "conta_bases.pl -i '{}' |grep '{}'".format(fasta_file, filtering)
            
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()
            if error or not output:
                return None
            if wanted=="contig_id":
              log.write("'{}' Find contig_id by size:{}\n".format(contig_id,cmd))
              new_contig_id = output.decode("utf-8").split(">")[-1].split("contig")[0].split("\n")[0].split("\t")[0]
              if "chromosome" in contig_id:
                new_contig_id=contig_id.replace("chromosome","")
              return new_contig_id
            elif wanted=="size":
              log.write("'{}' Find size by contig_id:{}\n".format(contig_id,cmd))
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
        log.write("'{}' Find contig in file: {}\n".format(contig_id,command))
    
        process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        output, err = process.communicate()
    
        if output:
            return True
        else:
            return False
    
    def count_sequences(filename):
        try:
            count = int(subprocess.check_output("grep -c '^>' " + filename, shell=True))
            return count
        except subprocess.CalledProcessError:
            return None
    
    def selectseq(organism_id, fasta_folder, contig_id,element,old):
        try:
            fasta_file = os.path.join(fasta_folder, organism_id + '.fasta')
            if element==True:
              output_file = os.path.join(fasta_folder, contig_id + '.fasta')
              cmd = "selectseq -l {} -s '{}' -o '{}'".format(contig_id, fasta_file, output_file)
              log.write("'{}' Extract contig: {}\n".format(old,cmd))
            else:
              output_file = os.path.join(fasta_folder, "{}-{}.fasta".format(organism_id,contig_id))
              cmd = "selectseq -l {} -c '{}' -o '{}'".format(contig_id,fasta_file,output_file)
              log.write("'{}' Subtract contig: {}\n".format(old,cmd))
            subprocess.check_call(cmd, shell=True)
            return True
        except subprocess.CalledProcessError as e:
            log.write("ERROR: " + str(e)+"\n")
            return False
    
    def blastn_element(outfile,element_fasta, organism_multifasta,contig_id):
        try:
          blastn_cline = NcbiblastnCommandline(query=element_fasta, subject=organism_multifasta, outfmt=6, out=outfile)
          stdout, stderr = blastn_cline()
          if os.path.isfile(outfile):
              if os.path.getsize(outfile) > 0:
                  with open(outfile, 'r') as f:
                      hits = [hit for hit in f.readlines()]
                  return hits[0].split()[1]
              else:
                log.write("'{}' ERROR BLASTN: No hits found!\n".format(contig_id))
                return None
          else:
              log.write("'{}' ERROR BLASTN: Result file was not created!\n".format(contig_id))
              return None
        except Exception as e:
          log.write("'{}' ERROR BLASTN: {}\n".format(contig_id,e))
    #      exc_type, exc_value, exc_traceback = sys.exc_info()
    #      formatted_lines = traceback.format_exc().splitlines()
    #      for i, line in enumerate(formatted_lines):
    #         log.write("{}: {}\n".format(i, line))
          return None
    
    def get_sequences(contig_id, fasta_folder, old,new):
        new_names=open(os.path.join(args.o,"selected_new_names.tab"), 'r+')
        coordinates = []
        sequence = ''
        contigs.append(contig_id)
        organism_name=None
        try:
            fasta_folder=os.path.abspath(fasta_folder)
            log.write("\nContig_id: {}\n".format(contig_id))
            log.write("'{}' fasta folder: {}\n".format(contig_id,fasta_folder))
            if not find_contig_in_file(contig_id.split("_")[0], args.old,contig_id):
              log.write("'{}' ERROR: Contig_id not found in table (-old)!\n".format(contig_id))
              if not "ERROR: Contig_id not found in table (-old)" in error.keys():
                  error["ERROR: Contig_id not found in table (-old)"]=[]
              error["ERROR: Contig_id not found in table (-old)"].append(contig_id)
            else:
              organism_name=get_value(old, contig_id.split("_")[0],"Organism name", "contig id")
              if organism_name==None:
                log.write("ERROR '{}': organism_name was not found!\n".format(contig_id))
                if not "ERROR: organism_name was not found" in error.keys():
                  error["ERROR: organism_name was not found"]=[]
                error["ERROR: organism_name was not found"].append(contig_id)
              else:
                log.write("'{}' Organism_name: {}\n".format(contig_id,organism_name))
                if not find_contig_in_file(contig_id, args.old,contig_id):
                  log.write("'{}' ERROR: Contig_id not found in table (-old)!\n".format(contig_id))
                  if not "ERROR: Contig_id not found in table (-old)" in error.keys():
                      error["ERROR: Contig_id not found in table (-old)"]=[]
                  error["ERROR: Contig_id not found in table (-old)"].append(contig_id)
                else:
                  contig_size=get_value(old, contig_id,"contig size","contig id")
                  if contig_size==None:
                    log.write("ERROR '{}': contig_size was not found!\n".format(contig_id))
                    if not "ERROR: contig_size was not found" in error.keys():
                      error["ERROR: contig_size was not found"]=[]
                    error["ERROR: contig_size was not found"].append(contig_id)
                  else:
                    log.write("'{}' Contig_size: {}\n".format(contig_id,contig_size)) 
                    old_organism_id=get_value(old, contig_id,"Organism id", "contig id")
                    if not old_organism_id==None:
                      log.write("'{}' Organism_id: {}\n".format(contig_id,old_organism_id))
                  if not organism_name==None:
                    organism_id=get_value(new, organism_name,"genome id","genome name")
                  if not organism_id==None and not old_organism_id==None:
                    if not old_organism_id==organism_id:
                      log.write("'{}' ATTENTION: organism_id has changed!\n".format(contig_id))
                  if organism_id==None:
                    log.write("ERROR '{}': Organism_id was not found!\n".format(contig_id))
                    if not "ERROR: Organism_id was not found" in error.keys():
                      error["ERROR: Organism_id was not found"]=[]
                    error["ERROR: Organism_id was not found"].append(contig_id)
                  else:
                    log.write("'{}' New organism_id: {}\n".format(contig_id,organism_id))
                    if not os.path.exists(fasta_folder.replace("\ "," ")):
                      os.makedirs(fasta_folder.replace("\ "," "))
                    fasta_folder=fasta_folder.replace("\\ "," ")
                    os.chdir(fasta_folder)
                    if not download_sequences(organism_id,contig_id):  
                      if not "ERROR: Sequenced genomes of the organism have not been downloaded" in error.keys():
                        error["ERROR: Sequenced genomes of the organism have not been downloaded"]=[]
                      error["ERROR: Sequenced genomes of the organism have not been downloaded"].append(contig_id)
                      log.write("'{}' ERROR: Sequenced genomes of the organism have not been downloaded!\n".format(contig_id))
                    elif not os.path.exists(os.path.join(fasta_folder,"{}.fasta".format(organism_id))):
                      if not "ERROR: Fasta file of sequenced genomes of the organism was not created" in error.keys():
                        error["ERROR: Fasta file of sequenced genomes of the organism was not created"]=[]
                      error["ERROR: Fasta file of sequenced genomes of the organism was not created"].append(contig_id)
                      log.write("'{}' ERROR: Fasta file of sequenced genomes of the organism was not created!\n".format(contig_id))
                    elif not os.path.getsize(os.path.join(fasta_folder,"{}.fasta".format(organism_id))) > 0:
                      if not "ERROR: Fasta file of sequenced genomes of the organism is empty" in error.keys():
                        error["ERROR: Fasta file of sequenced genomes of the organism is empty"]=[]
                      error["ERROR: Fasta file of sequenced genomes of the organism is empty"].append(contig_id)
                      log.write("'{}' ERROR: Fasta file of sequenced genomes of the organism is empty!\n".format(contig_id))
                    else:
                      if not contig_size==None:
                        new_contig_id=count_contig_bases(os.path.join(fasta_folder,"{}.fasta".format(organism_id)), contig_size,contig_id,"contig_id")   
                        if new_contig_id==None:
                          if not "ATTENTION: New contig_id not found based on size" in error.keys():
                            error["ATTENTION: New contig_id not found based on size"]=[]
                          error["ATTENTION: New contig_id not found based on size"].append(contig_id)
                          log.write("'{}' ATTENTION: New contig_id not found based on size!\n".format(contig_id)) 
                        else:
                          if "chromosome" in  new_contig_id:
                            new_contig_id=new_contig_id.replace("chromosome","")
                          if not "{}\t{}\n".format(contig_id,new_contig_id) in new_names.read():
                            new_names.write("{}\t{}\n".format(contig_id,new_contig_id))
                            new_names.close()
                            new_names=open(os.path.join(args.o,"selected_new_names.tab"), 'r+')
                          if not contig_id==new_contig_id:
                            log.write("'{}' New contig_id based on size: {}\n".format(contig_id,new_contig_id)) 
                      filename=''
                      namesize=9999
                      positive=None
                      for name in element_files.keys():
                        if contig_id.lower() in name.lower():
                          if len(name)<namesize:
                            filename=name
                            namesize=len(name)
                      if filename=='':
                        log.write("'{}' ERROR: Element's fasta was not found!\n".format(contig_id))
                        if not "ERROR: Element's fasta was not found" in error.keys():
                          error["ERROR: Element's fasta was not found"]=[]
                        error["ERROR: Element's fasta was not found"].append(contig_id)
                      else:
                        if not os.path.exists(os.path.join("/".join(fasta_folder.split("/")[0:-1]),'BLASTN').replace("\ "," ")):
                          os.mkdir(os.path.join("/".join(fasta_folder.split("/")[0:-1]),'BLASTN').replace("\ "," "))
                        outfile=os.path.join("/".join(fasta_folder.split("/")[0:-1]),'BLASTN','blastn_{}.txt'.format(contig_id))
                        positive=blastn_element(outfile, element_files[filename], os.path.join(fasta_folder,"{}.fasta".format(organism_id)),contig_id)
                        if positive==None:
                          if not "ERROR: No hits found" in error.keys():
                            error["ERROR: No hits found"]=[]
                          error["ERROR: No hits found"].append(contig_id)
                        else:
                          if "chromosome" in  positive:
                            positive=positive.replace("chromosome","")
                          if not "{}\t{}\n".format(contig_id,positive) in new_names.read():
                            new_names.write("{}\t{}\n".format(contig_id,positive))
                          log.write("'{}' New contig_id based on BLASTN: {}\n".format(contig_id,positive))
                          if new_contig_id==None or not new_contig_id==positive:
                              new_size=count_contig_bases(os.path.join(fasta_folder,"{}.fasta".format(organism_id)), positive,contig_id,"size")
                              if not new_size==None:
                                log.write("'{}' New contig size based on BLASTN: {}\n".format(contig_id,new_size))
                        if not positive==None and not positive==new_contig_id and not new_contig_id==None:
                          log.write("'{}' ERROR: Contig_id found by size and by BLASTN search was not a match!\n".format(contig_id))
                          if not "ERROR: Contig_id found by size and by BLASTN search was not a match" in error.keys():
                            error["ERROR: Contig_id found by size and by BLASTN search was not a match"]=[]
                          error["ERROR: Contig_id found by size and by BLASTN search was not a match"].append(contig_id)
                      
                      for cid in list(dict.fromkeys([positive,new_contig_id])):
                        if all(item is None for item in list(dict.fromkeys([positive,new_contig_id]))):
                          if contig_id not in not_created:
                            not_created.append(contig_id)
                          if contig_id in created:
                            created.remove(contig_id)
                        if not cid==None:
                          if "chromosome" in cid:
                            cid=cid.replace("chromosome","")
                          fasta_file = os.path.join(fasta_folder, organism_id + '.fasta')
                          if not find_contig_in_file(cid, fasta_file,contig_id):
                            if not "ERROR: Fasta file of sequenced organism genomes does not contain contig_id" in error.keys():
                              error["ERROR: Fasta file of sequenced organism genomes does not contain contig_id"]=[]
                            error["ERROR: Fasta file of sequenced organism genomes does not contain contig_id"].append(contig_id)
                            if contig_id not in not_created:
                              not_created.append(contig_id)
                            log.write("'{}' ERROR: Fasta file of sequenced organism genomes does not contain contig_id!\n".format(contig_id))
                          else:
                            if not selectseq(organism_id, fasta_folder, cid,True,contig_id):
                              if not "ERROR: selectseq exception when extracting contig" in error.keys():
                                error["ERROR: selectseq exception when extracting contig"]=[]
                              error["ERROR: selectseq exception when extracting contig"].append(contig_id)
                              log.write("'{}' ERROR: selectseq exception when extracting contig!\n".format(contig_id))  
                              if contig_id in created:
                                  created.remove(contig_id)
                              if contig_id not in not_created:
                                not_created.append(contig_id)
                            elif not os.path.exists(os.path.join(fasta_folder, cid + '.fasta')):
                              log.write("'{}' ERROR: Fasta file of the extracted contig was not created!\n".format(contig_id))  
                              if not "ERROR: Fasta file of the extracted contig was not created" in error.keys():
                                error["ERROR: Fasta file of the extracted contig was not created"]=[]
                              error["ERROR: Fasta file of the extracted contig was not created"].append(contig_id)
                              if contig_id in created:
                                  created.remove(contig_id)
                              if contig_id not in not_created:
                                not_created.append(contig_id)
                            elif not os.path.getsize(os.path.join(fasta_folder, cid + '.fasta')) > 0:
                              log.write("'{}' ERROR: Fasta file of the extracted contig is empty!\n".format(contig_id))  
                              if not "ERROR: Fasta file of the extracted contig is empty" in error.keys():
                                error["ERROR: Fasta file of the extracted contig is empty"]=[]
                              error["ERROR: Fasta file of the extracted contig is empty"].append(contig_id)
                              if contig_id in created:
                                  created.remove(contig_id)
                              if contig_id not in not_created:
                                not_created.append(contig_id)
                            else:
                              created.append(contig_id)
                          nseqs=count_sequences(fasta_file)
                          if nseqs==None:
                            if not "ERROR: Couldn't count the contigs in the organism's multifasta" in error.keys():
                              error["ERROR: Couldn't count the contigs in the organism's multifasta"]=[]
                            error["ERROR: Couldn't count the contigs in the Fasta file of the organism's multifasta"].append(contig_id)
                            log.write("'{}' ERROR: Couldn't count the contigs in the Fasta file of the organism's multifasta!\n".format(contig_id))
                          elif nseqs==0:
                            if not "ERROR: The organism's multifasta has no contig" in error.keys():
                              error["ERROR: The organism's multifasta has no contig"]=[]
                            error["ERROR: The organism's multifasta has no contig"].append(contig_id)
                            log.write("'{}' ERROR: The organism's multifasta has no contig!\n".format(contig_id))
                            if contig_id in created:
                              created.remove(contig_id)
                            if contig_id not in not_created:
                              not_created.append(contig_id)
                          elif nseqs==1:
                            if not "ATTENTION: The organism's multifasta has only one contig" in error.keys():
                              error["ATTENTION: The organism's multifasta has only one contig"]=[]
                            error["ATTENTION: The organism's multifasta has only one contig"].append(contig_id)
                            log.write("'{}' ATTENTION: The organism's multifasta has only one contig!\n".format(contig_id))
                          else:
                            if not selectseq(organism_id, fasta_folder, cid,False,contig_id):
                              if not "ERROR: selectseq exception when subtracting contig" in error.keys():
                                error["ERROR: selectseq exception when subtracting contig"]=[]
                              error["ERROR: selectseq exception when subtracting contig"].append(contig_id)
                              log.write("'{}' ERROR: selectseq exception when subtracting contig!\n".format(contig_id))
                              if contig_id not in not_created:
                                  not_created.append(contig_id)
                              if contig_id in created:
                                  created.remove(contig_id)
                            elif not os.path.exists(os.path.join(fasta_folder, "{}-{}.fasta".format(organism_id,cid))):
                              log.write("'{}' ERROR: Fasta file of the subtracted contig was not created!\n".format(contig_id))  
                              if not "ERROR: Fasta file of the subtracted contig was not created" in error.keys():
                                error["ERROR: Fasta file of the subtracted contig was not created"]=[]
                              error["ERROR: Fasta file of the subtracted contig was not created"].append(contig_id)
                              if contig_id not in not_created:
                                  not_created.append(contig_id)
                              if contig_id in created:
                                  created.remove(contig_id)
                            elif not os.path.getsize(os.path.join(fasta_folder, "{}-{}.fasta".format(organism_id,cid))) > 0:
                              log.write("'{}' ERROR: Fasta file of the subtracted contig is empty!\n".format(contig_id))  
                              if not "ERROR: Fasta file of the subtracted contig is empty" in error.keys():
                                error["ERROR: Fasta file of the subtracted contig is empty"]=[]
                              error["ERROR: Fasta file of the subtracted contig is empty"].append(contig_id)
                              if contig_id not in not_created:
                                  not_created.append(contig_id)
                              if contig_id in created:
                                  created.remove(contig_id)
                            else:
                              if not contig_id in created:
                                  created.append(contig_id)
                      os.chdir(call)
                      
    
        except Exception as e:
            if "gres Exception" not in error.keys():
              error["gres Exception"]=[]
            error["gres Exception"].append(contig_id)
            print("ERROR '{}': {}".format(contig_id, e))
            exc_type, exc_value, exc_traceback = sys.exc_info()
            formatted_lines = traceback.format_exc().splitlines()
            for i, line in enumerate(formatted_lines):
                print("{}: {}".format(i, line))
        new_names.close()
    
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
        new_names.write("Old\tNew\n")
      except Exception as e:
        print("ERROR: Selected_new_names.tab file was not written!")
        exit()
    new_names.close()
    
    log=open(os.path.join(args.o,"file.log"), 'a')
    for family in sorted(contig_list.keys()):
        log.write("\nFAMILY: {}\n".format(family.upper()))
        os.chdir(call)
        if not family=="":
          fasta_folder=os.path.join(args.o,family)
        else:
          fasta_folder=args.o
        fasta_folder=os.path.join(fasta_folder,"FASTA")
        if not os.path.exists(fasta_folder):
          os.makedirs(fasta_folder.replace("\ "," "))
        contig_list[family]=list(set(contig_list[family]))
        for contig_id in sorted(contig_list[family]):
          get_sequences(contig_id,fasta_folder,old,new)
          log.close()
          log=open(os.path.join(args.o,"file.log"), 'a')
except Exception as e:
    if "gres Exception" not in error.keys():
      error["gres Exception"]=[]
    error["gres Exception"].append(contig_id)
    print("ERROR '{}': {}".format(contig_id, e))
    exc_type, exc_value, exc_traceback = sys.exc_info()
    formatted_lines = traceback.format_exc().splitlines()
    for i, line in enumerate(formatted_lines):
        print("{}: {}".format(i, line))
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
contigs=list(set(contigs))
if not len(contigs)==0:
    log.write("\n{} contigs were processed: {}\n".format(len(contigs),sorted(contigs)))
    print("{} contigs were processed!".format(len(contigs)))
print("End")
log.write("\nEnd\n")
log.close()
new_names.close()