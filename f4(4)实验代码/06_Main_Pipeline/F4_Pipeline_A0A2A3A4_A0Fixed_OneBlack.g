# 入口脚本用途:
# - 定义 a0 固定黑，且 a2/a3/a4 三选一黑 的实验
# - target_color_patterns 中每一项对应一类允许配置
# - Main_Pipeline.g 会遍历并筛出满足模式的子系统
GlobalOutFile := "../f4(4)_a0a2a3a4_a0固定黑且a2a3a4任一单黑_实验运行结果.txt";;
GlobalIntroText := "研究对象: F4(4) 实验 -> 选取[a0,a2,a3,a4]，其中a0固定涂黑，且a2/a3/a4中任一位置涂黑";;

GlobalExperiments := [
    rec(
        desc := "F4(4) 实验: 选取[a0,a2,a3,a4]，a0固定涂黑，且a2/a3/a4中任意一个涂黑",
        roots := rec(
            mode := "SCAN",
            selected_pos_labels := ["a0","a2","a3","a4"],
            target_color_patterns := [
                [1,1,0,0],
                [1,0,1,0],
                [1,0,0,1]
            ]
        ),
        colors := fail
    )
];;

Read("Main_Pipeline.g");;
