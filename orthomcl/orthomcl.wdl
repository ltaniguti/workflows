task compliantFasta {
    File proteome
    String identifier

    command {
        /software/orthomclSoftware-v2.0.9/bin/orthomclAdjustFasta ${identifier} ${proteome} 1
    }

    runtime {
        docker: "gcr.io/esalq-carvao/orthomcl:v1"
        memory: "1 GB"
        cpu: "1"
        disks: "local-disk " + 10 + " HDD"
        preemptible: 3  
    }

    output {
        File fasta = "${identifier}.fasta"
    }
}

task PrepareOrthoMCL {
  Array[File] compliant_fastas

  command {
    echo [mysqld] >> /etc/alternatives/my.cnf
    echo myisam_sort_buffer_size=3G >> /etc/alternatives/my.cnf
    echo myisam_max_sort_file_size=15000 >> /etc/alternatives/my.cnf
    echo [mysqldump] >> /etc/alternatives/my.cnf
    /etc/init.d/mysql start

    mysql -u root --execute "create database orthomcl;"

    echo -e "dbVendor=mysql
dbConnectString=dbi:mysql:orthomcl:localhost:3307
dbLogin=root
dbPassword=NONE
similarSequencesTable=SimilarSequences
orthologTable=Ortholog
inParalogTable=InParalog
coOrthologTable=CoOrtholog
interTaxonMatchView=InterTaxonMatch
percentMatchCutoff=50
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE" > orthomcl.config

    /software/orthomclSoftware-v2.0.9/bin/orthomclInstallSchema orthomcl.config


    mkdir compliantFasta
    cp ${sep=" " compliant_fastas} compliantFasta/

    /software/orthomclSoftware-v2.0.9/bin/orthomclFilterFasta compliantFasta/ 10 20
    makeblastdb -in goodProteins.fasta -dbtype prot -out all_prot.db
    blastp -seg yes -evalue 1e-5 -db all_prot.db -query goodProteins.fasta -outfmt 6 -out all-vs-all.tsv -num_threads 8
    /software/orthomclSoftware-v2.0.9/bin/orthomclBlastParser all-vs-all.tsv compliantFasta >> similarSequences.txt
    /software/orthomclSoftware-v2.0.9/bin/orthomclLoadBlast orthomcl.config similarSequences.txt
    /software/orthomclSoftware-v2.0.9/bin/orthomclPairs orthomcl.config orthomcl_pairs.log cleanup=no
    /software/orthomclSoftware-v2.0.9/bin/orthomclDumpPairsFiles orthomcl.config
    mcl mclInput --abc -I 1.5 -o mclOutput
    /software/orthomclSoftware-v2.0.9/bin/orthomclMclToGroups gus 1 < mclOutput > groups.txt
  }

  runtime {
    docker: "gcr.io/esalq-carvao/orthomcl:v1"
    memory: "8 GB"
    cpu: "4"
    disks: "local-disk " + 10 + " HDD"
    preemptible: 3  
  }

  output {
    File groups = "groups.txt"
    File blast = "all-vs-all.tsv"
  }
}


workflow Orthomcl {
    Array[File] proteomes
    Array[String] identifiers

    Array[Pair[File, File]] pairs = zip(proteomes, identifiers)

    scatter (pair in pairs) {
        call compliantFasta {
            input:
                proteome=pair.left,
                identifier=pair.right
        }
    }

    call PrepareOrthoMCL {
        input:
            compliant_fastas=compliantFasta.fasta
    }

    output {
        File final_groups = PrepareOrthoMCL.groups
        File blast_result = PrepareOrthoMCL.blast
    }
}
