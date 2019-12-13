# Almost all eukaryotic secreted proteins contain a signal peptide at the N-terminus 
# that directs proteins to the rough ER and the Golgi complex (Blobel and Dobberstein, 
# 1975; von Heijne, 1990). The signal peptide, typically 15 – 30 amino acids long, is 
# cleaved off during translocation across the membrane. While some proteins without an 
# N-terminal signal peptide can be found in the ER and the Golgi, over 90% of human 
# secreted proteins contain classical N-terminal signal peptides (Scott et al., 2004)


# For SignalP prediction, only entries that are predicted to have a "mostly likely 
# cleavage site" by SignalP-NN algorithm and a &signal peptide" by SignalP-HMM algorithm 
# are considered to be true "positives" (Bendtsen et al., 2004) using the N-terminal 70 amino 
# acids
task SignalPeptide {
  File proteins

  command {
    signalp -fasta ${proteins} -org euk -prefix out
  }

  runtime {
    docker: "gcr.io/taniguti-backups/bio/signalp:latest"
    memory: "1 GB"
    cpu: "1"
    disks: "local-disk " + 10 + " HDD"
    preemptible: 3
  }

  output {
    File signal = "out_summary.signalp5"
  }

}


# For predicting membrane proteins using TMHMM, the entries having membrane 
# domains not located within the N-terminus (the first 70 amino acids) were 
# treated as real membrane proteins since our analysis showed treating entries 
# having a single transmembrane domain located in the N-terminus significantly 
# decreased the prediction accuracy.
#
# allowing one transmembrane (TM) domain in the first 60 aa
task TransmembraneDomain {
  File proteins

  command {
    cat ${proteins} | tmhmm > transmembrane.txt
  }

  runtime {
    docker: "taniguti/tmhmm:latest"
    memory: "1 GB"
    cpu: "1"
    disks: "local-disk " + 10 + " HDD"
    preemptible: 3
  }

  output {
    File transmembrane = "transmembrane.txt"
  }

}

# removing proteins with ER targeting sequence (Prosite: PS00014)
task PrositeScan {
  File proteins

  command {
    perl /ps_scan/ps_scan.pl -d /ps_scan/prosite.dat -p PS00014 ${proteins} > ER.txt
  }

  runtime {
    docker: "taniguti/psscan:latest"
    memory: "1 GB"
    cpu: "1"
    disks: "local-disk " + 10 + " HDD"
    preemptible: 3
  }

  output {
    File er = "ER.txt"
  }
}


task selectSecretome {
  File signalp
  File transmembrane
  File endoplasmatic_reticulum
  File proteome

  command <<<
    python3.5 <<CODE

    import sys
    from Bio import SeqIO


    def read_file(path):
        lines = []
        with open(path) as f:
            for l in f:
                if l.startswith("#"):
                    continue
                lines.append(l.strip().split())
        return lines


    def parse_signalp(items):
        if items[1] != "OTHER":
            ident = items[0]
            signal = "Y"
            return ident, signal
        return items[0], "N"

    def parse_tmhmm(items):
        prot, __, tipo, start, end = items
        return prot, tipo, start, end

    def interpret_tmhmm(regions):
        n_terminal = True
        for cat, start, end in regions:
            if cat != 'outside':
                if end < 60:
                    continue
                n_terminal = False
        return n_terminal

    def directed_to_er(path):
        idents = []
        with open(path) as f:
            for l in f:
                if l.startswith(">"):
                    idents.append(l.strip().split(" ")[0].replace(">", ""))
        return idents

    def select_secretome_subset(wanted, whole_proteome):
        seqiter = SeqIO.parse(whole_proteome, 'fasta')                                    
        SeqIO.write((seq for seq in seqiter if seq.id in wanted), sys.stdout, "fasta")

    file_signalp = "${signalp}"
    file_tmhmm = "${transmembrane}"
    file_endoplasmatic_reticulum = "${endoplasmatic_reticulum}"
    proteome = "${proteome}"
    

    er_pred = directed_to_er(file_endoplasmatic_reticulum)
    signalp_lines = read_file(file_signalp)
    signal_pred = {}
    for l in signalp_lines:
        prot, pred = parse_signalp(l)
        signal_pred[prot] = pred

    tmhmm_lines = read_file(file_tmhmm)
    tmhmm_pred = {}
    for l in tmhmm_lines:
        prot, tipo, start, end = parse_tmhmm(l)
        tmhmm_pred.setdefault(prot, []).append((tipo, int(start), int(end)))

    fp_tmhmm = []
    for i, v in tmhmm_pred.items():
        is_n_terminal_only = interpret_tmhmm(v)
        if is_n_terminal_only:
            fp_tmhmm.append(i)
    
    selected = []
    for k, v in signal_pred.items():
        if v == "Y" and k in fp_tmhmm and k not in er_pred:
            selected.append(k)


    select_secretome_subset(selected, proteome)

    CODE
  >>>

  runtime {
    docker: "biocontainers/biopython:v1.68dfsg-3-deb-py3_cv1"
    memory: "1 GB"
    cpu: "1"
    disks: "local-disk " + 10 + " HDD"
    preemptible: 3
  }

  output {
    File secretome = stdout()
  }
}

task localize {
  File proteins

  command {
    python /LOCALIZER_1.0.4/Scripts/LOCALIZER.py -e -i ${proteins} > location.txt
  }

  runtime {
    docker: "taniguti/localizer:latest"
    memory: "1 GB"
    cpu: "1"
    disks: "local-disk " + 10 + " HDD"
    preemptible: 3
  }

  output {
    File location = "location.txt"
  }
}

task EffectorP {
  File secretome

  command {
    python /EffectorP_2.0/Scripts/EffectorP.py -s -E effectors.fa -i "${secretome}"
  }

  runtime {
    docker: "taniguti/effectorp2:latest"
    memory: "1 GB"
    cpu: "1"
    disks: "local-disk " + 10 + " HDD"
    preemptible: 3
  }

  output {
    File effectors = stdout()
    File fasta = "effectors.fa"
  }
}

task RunPhibase {
  File phibase
  File secretome

  command {
    makeblastdb -dbtype prot -in ${phibase}
    blastp -max_target_seqs 5 -evalue 1e-5 -query ${secretome}  \
           -db ${phibase} \
           -outfmt "6 qseqid sseqid qstart qend sstart send evalue qcovs" -out resultados.txt
  }

  runtime {
    docker: "taniguti/orthomcl:latest"
    memory: "8 GB"
    cpu: "4"
    disks: "local-disk " + 10 + " HDD"
    preemptible: 3  
  }

  output {
    File blast_result = "resultados.txt"
  }
}

task RunCazy {
  File secretome
  File cazydb

  command {
    hmmsearch --tblout hmmer_out.txt ${cazydb} ${secretome}
  }

  runtime {
    docker: "dockerbiotools/hmmer:latest"
    memory: "8 GB"
    cpu: "4"
    disks: "local-disk " + 10 + " HDD"
    preemptible: 3   
  }

  output {
    File hmmer_result = "hmmer_out.txt"
  }
}

workflow Secretome {
  File proteins
  File phibase
  File cazydb

  call SignalPeptide {
    input: proteins=proteins
  }

  call TransmembraneDomain {
    input: proteins=proteins
  }

  call PrositeScan {
    input: proteins=proteins
  }

  call selectSecretome {
    input:
      signalp=SignalPeptide.signal,
      transmembrane=TransmembraneDomain.transmembrane,
      endoplasmatic_reticulum=PrositeScan.er,
      proteome=proteins
  }
  
  call localize {
    input: proteins=selectSecretome.secretome
  }

  call EffectorP {
    input:
      secretome=selectSecretome.secretome
  }

  call RunPhibase {
    input:
      phibase=phibase,
      secretome=secretome
  }

  call RunCazy {
    input:
      secretome=secretome,
      cazydb=cazydb
  }

  output {
    File sig = SignalPeptide.signal
    File trans = TransmembraneDomain.transmembrane
    File pscan = PrositeScan.er
    File local = localize.location
    File secretome = selectSecretome.secretome
    File effectors = EffectorP.effectors
    File blast_phibase = RunPhibase.blast_result
    File cazy_result = RunCazy.hmmer_result
  }
}
