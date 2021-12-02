# Miniconda usage

## 1 Download

```
wget -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda2-4.5.11-Linux-x86_64.sh
```

## 2 Install

```bash
bash Miniconda2-4.5.11-Linux-x86_64.sh  # uname -a
# 按enter--三下空格--输入yes--按enter--输入yes
source ~/.bashrc # 激活配置
```

## 3 Configuration

do only once

```
conda config --append channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free \
conda config --append channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge \
conda config --append channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda \
conda config --set show_channel_urls yes

conda config --append default_channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
conda config --append default_channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda config --append default_channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
```

执行完上述命令后，会生成\~/.condarc(这里面包含了用户对miniconda的配置信息)

```
vim ~/.condarc
```

check config info

```
conda config --show
```

## 4 Create environment

```
# 创建名为rna的软件安装环境，同时安装python=2版本的软件
conda create -n rna python=2
# 出现三个done
```

```
# 查看当前conda环境
# 可以看到成功建立的rna
conda info --envs
```

```
# 激活/进入conda的rna环境，避免安装软件时安装到大环境
source activate rna
# 小环境创建成功，可以随便安装软件到小环境里啦
```

## 5 Install package

```
# 安装 sra-tools软件
conda search sra-tools  
conda install -y sra-tools  # done正确安装，且能调出软件help
......
source deactivate # 退出当前环境
```

## 6 Others

```
# 查看环境名：
conda info --envs 或conda info -e
# 查看已安装软件列表：
conda list
# 退出环境：
source deactivate
# 更新：
conda update python 
# (conda将python等软件都视为package)
# 假设当前环境是python 3.4, conda会将python升级为3.4.x系列的当前最新版本

# 删除全部packages
conda remove --name/-n wes --all
# 删除某个packages
source activate wes
conda remove multiqc
# 或直接指明name
conda remove -n wes numpy
```
