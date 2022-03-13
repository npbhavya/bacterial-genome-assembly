    #!/bin/bash

    #SBATCH --job-name=guppy_basecaller
    #SBATCH --mail-user=email
    #SBATCH --mail-type=ALL
    #SBATCH --output=%x-%j.out.txt
    #SBATCH --error=%x-%j.err.txt
    #SBATCH --time=1-0
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50G
    #SBATCH --gres=gpu:1
    #SBATCH --partition=hpc_gpu

    #Runs guppy basecaller GPU version on a GPU node
    module load cuda10.0/nsight/10.0.130
    export PATH=~/ont-guppy/bin:$PATH

    cd <path to run the script and find the files>
    guppy_basecaller --input_path <directory with fast5 files> --save_path <output path> --flowcell FLO-FLG001 --kit SQK-RBK004 --device auto
