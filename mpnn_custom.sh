#!/bin/bash -l
#SBATCH -J MPNN_for_2spots
#SBATCH -n 3
#SBATCH -N 1
#SBATCH -p helium
#SBATCH --mem=32G
### # SBATCH --gres=gpu:1
#SBATCH -o 2mpnn1.log
#SBATCH -e 2mpnn1.log
#SBATCH --qos=qcpu

# 激活环境
conda activate /xcfhome/yzmeng/miniconda3/envs/dl_binder_design

# 检查参数数量
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

# 读取参数
pdbdir=$1
outpdbdir=$2

# 如果目录不以 / 结尾，就补上
[[ "${pdbdir}" != */ ]] && pdbdir="${pdbdir}/"
[[ "${outpdbdir}" != */ ]] && outpdbdir="${outpdbdir}/"

# 若目录不存在则创建
if [ ! -d "$pdbdir" ]; then
    echo "Input directory does not exist, creating: $pdbdir"
    mkdir -p "$pdbdir"
fi

if [ ! -d "$outpdbdir" ]; then
    echo "Output directory does not exist, creating: $outpdbdir"
    mkdir -p "$outpdbdir"
fi

# 打印最终路径，方便调试
echo "Using input directory:  $pdbdir"
echo "Using output directory: $outpdbdir"

# 运行主程序
python /xcfhome/yzmeng/tools/design/dl_binder_design_3/mpnn_fr/dl_interface_design.py \
    -pdbdir "$pdbdir" \
    -outpdbdir "$outpdbdir" \
    -relax_cycles 0 \
    -seqs_per_struct 50
