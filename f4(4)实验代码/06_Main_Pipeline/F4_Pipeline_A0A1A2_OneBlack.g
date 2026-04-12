# 入口脚本用途:
# - 定义单个 F4 实验场景（a0,a1,a2 三点中恰一黑）
# - 通过 GlobalOutFile 指定结果写入目标
# - 通过 GlobalExperiments 传参给 Main_Pipeline.g
# - 本文件不含算法逻辑，仅负责“配置注入 + 调主流程”
GlobalOutFile := "../f4(4)_a0a1a2_任一单黑_实验运行结果.txt";;
GlobalIntroText := "研究对象: F4(4) 实验 -> 选取[a0,a1,a2] 且恰有一个黑点";;

GlobalExperiments := [
    rec(
        desc := "F4(4) 实验: 选取[a0,a1,a2]，三个根中任意一个涂黑",
        roots := rec(
            mode := "SCAN",
            selected_pos_labels := ["a0","a1","a2"],
            required_black_count := 1
        ),
        colors := fail
    )
];;

Read("Main_Pipeline.g");;
