# 入口脚本用途:
# - 定义 a0 固定黑，且 a2/a3 二选一黑 的实验
# - 通过 target_color_patterns 明确列出允许的两种涂黑向量
# - 其余流程（模块2~5）均复用 Main_Pipeline.g
GlobalOutFile := "../f4(4)_a0a2a3_a0固定黑且a2a3任一单黑_实验运行结果.txt";;
GlobalIntroText := "研究对象: F4(4) 实验 -> 选取[a0,a2,a3]，其中a0固定涂黑，且a2/a3中任一位置涂黑";;

GlobalExperiments := [
    rec(
        desc := "F4(4) 实验: 选取[a0,a2,a3]，a0固定涂黑，且a2/a3中任意一个涂黑",
        roots := rec(
            mode := "SCAN",
            selected_pos_labels := ["a0","a2","a3"],
            target_color_patterns := [
                [1,1,0],
                [1,0,1]
            ]
        ),
        colors := fail
    )
];;

Read("Main_Pipeline.g");;
