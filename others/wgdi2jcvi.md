# wgdi2jcvi

使用wgdi结果绘制jcvi的两种图：`Macrosynteny getting fancy` 和`Microsynteny getting fancy`。

## （一）Macrosynteny 大尺度共线性

### 1.数据准备及基础运行方式

#### (1)`bed`文件

六列：染色体表示，起始，终止，id，占位0，方向。

![image-20250110092857938](C:\Users\yuxia\AppData\Roaming\Typora\typora-user-images\image-20250110092857938.png)

#### (2）`seqids`文件

包含需要绘制的染色体标识信息，每个物种一行，顺序与layout中的track顺序保持一致。

![image-20250110093234799](C:\Users\yuxia\AppData\Roaming\Typora\typora-user-images\image-20250110093234799.png)

#### (3）`layout`文件

绘图布局文件，包含两部分信息

第一个#内(控制染色体的排布)

y轴，x轴起始，x轴终止，旋转，颜色，文本标签，垂直对齐方式，bed文件。

第二个#内(控制染色体内的连线)

固定占位e，连线的track1，连线的track2，simple文件

![image-20250110093606839](C:\Users\yuxia\AppData\Roaming\Typora\typora-user-images\image-20250110093606839.png)

#### (4) 运行方式

```shell
# 两个位置参数：seqids文件 和 layout文件 的路径。
python -m jcvi.graphics.karyotype seqids layout
```

### 2. 数据生成方式

(1) 生成`bed`文件

```shell
python wgdi2jcvi.py gff2bed -gff input_gff_file -bed output_bed_file

# 
# input_gff_file         wgdi支持的gff格式
# output_bed_file        jcvi支持的bed格式

```

(2)生成`simple`文件

```shell
python wgdi2jcvi.py block2simple -c input_collinearity_file -o output_simple_file [--filter filter_blockid_file]

# -c input_collinearity_file      wgdi生成的collinearity文件
# -o output_simple_file           jcvi支持的simple文件
# --filter filter_blockid_file    需要保留的block，可使用wgdi生成的筛选后的blockinfo文件的ID列。

# 举例1 不使用block筛选
python wgdi2jcvi.py block2simple -c coca_schi.collinearity.txt -o coca.schi.simple

# 举例2 使用block筛选
# coca_chi_filter_blockid.txt 来源于最终版blockinfo的ID列。
python wgdi2jcvi.py block2simple -c coca_schi.collinearity.txt -o coca.schi.simple --filter coca_chi_filter_blockid.txt
```

## (二) Microsynteny 小尺度共线性

### 2.数据准备及基础运行方式

#### (1)`blocks`文件

包含共线的基因对，可以从wgdi的alignment中截取一部分。

![image-20250110102201648](C:\Users\yuxia\AppData\Roaming\Typora\typora-user-images\image-20250110102201648.png)

#### (2）`block.layout`文件

绘图布局文件，包含两部分信息

第一个#内(控制基因对的排布)，每个物种的基因排布占一行

```
# x,   y, rotation,   ha,     va,   color, ratio,            label
x轴, y轴, 旋转, 水平对齐方式， 垂直对齐方式， 缩放比例，  文本标签
```

第二个#内(控制基因对间的连线)

固定占位e，连线的track1，连线的track2

![image-20250110104009448](C:\Users\yuxia\AppData\Roaming\Typora\typora-user-images\image-20250110104009448.png)

#### (3) 合并的`bed`文件

```shell
cat xxx.bed yyy.bed > xxx_yyy.bed
```

#### (4) 运行方式

```shell
python -m jcvi.graphics.synteny blocks xxx_yyy.bed blocks.layout

# blocks         想要展示局部blocks文件
# xxx_yyy.bed    合并的bed文件
# blocks.layout  blocks.layout文件 

# 范例：
python -m jcvi.graphics.synteny blocks grape_peach.bed blocks.layout
```

