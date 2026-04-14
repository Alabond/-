#############################################################################
# 文件: G2_Pipeline.g
# 描述: 使用 Main_Pipeline.g 的通用流程，配置为运行 G2 的两个实验
# 作用: 作为 G2 专用入口文件，负责预设全局变量，然后调用通用主流程。
# 说明: 若把 Main_Pipeline.g 看成“引擎”，那么本文件就是“G2 的启动配置”。
#############################################################################

# 配置输出文件与简介文本
GlobalOutputDir := "..";;
GlobalOutFile := Concatenation(GlobalOutputDir, "/g2(2)实验运行结果.txt");;
GlobalIntroText := "研究对象: G2(2) 两组实验 -> A2(一黑一白), A1+A1(两黑)";;
GlobalAmbientType := "G2";;
GlobalRequireModule4AllPass := false;;

# 配置实验：
# - 实验 1: A2 (长-长, 一黑一白)
# - 实验 2: A1 + A1 (两黑，正交)
GlobalExperiments := [
    rec(
        desc := "G2 实验 1: 长-长 一黑一白 -> A2",
        roots := rec(
            mode := "SCAN",
            template := [[0,1], [0,1]],   # 两个长根
            selected_pos_labels := ["a0", "a2"],
            target_black := 1,
            target_white := 1
        ),
        colors := fail
    ),
    rec(
        desc := "G2 实验 2: 两黑 正交 -> A1 + A1",
        roots := rec(
            mode := "SCAN",
            template := [[1,0], [0,1]],   # 模板占位，不影响判定
            target_black := 2,
            target_white := 0
        ),
        colors := fail
    )
];;

# 调用通用主流程
# 所有 G2 专用参数都在上面先设好，然后统一交给 Main_Pipeline.g 执行。
Read("Main_Pipeline.g");;
