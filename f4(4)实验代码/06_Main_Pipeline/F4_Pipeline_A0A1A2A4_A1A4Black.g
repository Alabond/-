# 入口脚本用途:
# - 定义单个 F4 实验场景（a0,a1,a2,a4，且 a1/a4 固定为黑）
# - 通过 selected_pos_labels + target_black_labels 精确限定涂黑位置
# - 将配置交给 Main_Pipeline.g 执行模块1~5全流程
GlobalOutFile := "../f4(4)_a0a1a2a4_a1a4涂黑_实验运行结果.txt";;
GlobalIntroText := "研究对象: F4(4) 实验 -> 选取[a0,a1,a2,a4]，固定a1,a4涂黑";;

GlobalExperiments := [
    rec(
        desc := "F4(4) 实验: 选取[a0,a1,a2,a4] 且 a1,a4 涂黑",
        roots := rec(
            mode := "SCAN",
            selected_pos_labels := ["a0","a1","a2","a4"],
            target_black_labels := ["a1","a4"]
        ),
        colors := fail
    )
];;

Read("Main_Pipeline.g");;
