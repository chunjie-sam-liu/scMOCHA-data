import "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/scMOCHA.wdl" as sub

workflow scMOCHABatch {

  String version = "v0.2.1"

  # Cell ranger inputs
  File output_id_list
  Array[String] output_ids = read_lines(output_id_list)

  File fastqs_list
  Array[String] fastqss = read_lines(fastqs_list)

  File sample_id_list
  Array[String] sample_ids = read_lines(sample_id_list)

  String chemistry = "auto"

  String transcriptome = "/home/liuc9/data/refdata/mgatk_index/Human"
  File rCRS = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/fasta/rCRS.MT.fasta"
  File mt_exons_df = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/fasta/mt_exons.df.rds.gz"
  File mt_features_gmoviz = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/fasta/mt_features.grange.gmoviz.rds.gz"

  File output_dir_list
  Array[String] output_dirs = read_lines(output_dir_list)

  # mgatk inputs
  String chrM = "MT"
  Int low_coverage_threshold = 10

  # cell_cluster_annotation inputs
  Int npcs = 10
  Float reso = 0.1
  String cellrefname
  String celllevel
  Int nFeature_RNA_min = 500
  Int nFeature_RNA_max = 6000
  Float percent_mt_max = 75
  Float percent_ribo_max = 50
  Float percent_Lagest_Gene_max = 50
  String x10_version = "v3"

  # Runtime attributes
  String memory = "50 GB"
  Int boot_disk_size_gb = 12
  String disk_space = "50"
  Int cpu = 10
  Boolean use_ssd = false

  # docker image
  String scmocha_version = "latest"
  String docker = "chunjiesamliu/scmocha"
  String partition = "defq"
  String account = "liuc9"
  File IMAGE = "/scr1/users/liuc9/sif/scmocha_latest.sif"

  File perlscript = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/bin/get_variants_info.pl"
  File jar_path = "/scr1/users/liuc9/tools/haplogrep3"  # /opt/haplogrep3/haplogrep3.jar
  File sqlite_path = "/mnt/isilon/xing_lab/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3"

  String bindir = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/bin"
  String conda_root = "/home/liuc9/tools/anaconda3"
  String conda_env = "scmocha"

  scatter (idx in range(length(output_ids))) {
    call sub.scMOCHA {
      input:
        version = version,

        # Cell ranger inputs
        output_id = output_ids[idx],
        fastqs = fastqss[idx],
        sample_id = sample_ids[idx],
        chemistry = chemistry,

        transcriptome = transcriptome,
        rCRS = rCRS,
        mt_exons_df = mt_exons_df,
        mt_features_gmoviz = mt_features_gmoviz,

        output_dir = output_dirs[idx],

        # mgatk inputs
        chrM = chrM,
        low_coverage_threshold = low_coverage_threshold,

        # cell_cluster_annotation inputs
        npcs = npcs,
        reso = reso,
        cellrefname = cellrefname,
        celllevel = celllevel,
        nFeature_RNA_min = nFeature_RNA_min,
        nFeature_RNA_max = nFeature_RNA_max,
        percent_mt_max = percent_mt_max,
        percent_ribo_max = percent_ribo_max,
        percent_Lagest_Gene_max = percent_Lagest_Gene_max,
        x10_version = x10_version,

        # Runtime attributes
        memory = memory,
        boot_disk_size_gb = boot_disk_size_gb,
        disk_space = disk_space,
        cpu = cpu,
        use_ssd = use_ssd,

        # docker image
        scmocha_version = scmocha_version,
        docker = docker,
        partition = partition,
        account = account,
        IMAGE = IMAGE,

        perlscript = perlscript,
        jar_path = jar_path,
        sqlite_path = sqlite_path,

        bindir = bindir,
        conda_root = conda_root,
        conda_env = conda_env
    }
  }


  output {
    Array[File] output_dir_tar_gzs = scMOCHA.output_dir_tar_gz
  }

}
