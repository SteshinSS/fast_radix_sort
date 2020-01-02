I have written boost:spreadsort-like sort, but it works 40% faster. The russian post is about how i did it.

# Сделаем spreadsort быстрым снова
## Задача
`Boost::spreadsort` -- это гибридный алгоритм in-place сортировки. Он запуксает radix/bucket-sort и переходит на `std::sort`, когда массив мал. Я реализовал ту же идею, но улучшил производительность. 

Всегда есть соблазн написать быстрый код, сузив исходную задачу, но мне это не подходит. Я постарался сохранить универсальность кода и, по-моему, немного увеличил ее. Формально, я пишу сортировку, которая, как и `spreadsort`:
- in-place
- non-stable
- single-thread
- реализует идею radix-sort

Я прочитал код `boost::spreadsort` и заметил еще несколько ограничений. Буст использует старшие биты ключа в качестве номера бакета. Это значит, что алгоритм расчитывает на равномерное распределение ключей. Именно на таких ключах я повторил результат из вики: двухкратное ускорение, по сравнению с `std::sort`. И еще: чтобы ускорение вышло именно за счет поразрядной сортировки, мы не будем использовать никакую сортировку, кроме `std::sort`.
- равномерное распределение ключей
- можно использовать только `std::sort`

Всегда можно выиграть, затачиваясь под свое железо. Для честного сравнения я сначала вручную поигрался с параметрами `spreadsort`, а потом запустил оптимизированный монте-карло поиск самых эффективных параметров. Через полтора часа и 20'000 итераций, я нашел несколько хороших наборов, но они почти не отличались от стандартных, поэтому я ничего не менял.
