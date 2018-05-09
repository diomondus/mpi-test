## mpi-test

#Настройки для macOS:

`brew install openmpi`

#MPI:

MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
IN buf	-	адрес начала расположения пересылаемых данных;
IN count	-	число пересылаемых элементов;
IN datatype	-	тип посылаемых элементов;
IN dest	-	номер процесса-получателя в группе, связанной с коммуникатором comm;
IN tag	-	идентификатор сообщения (аналог типа сообщения функций nread и nwrite PSE nCUBE2);
IN comm	-	коммуникатор области связи.


int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
OUT	buf	-	адрес начала расположения принимаемого сообщения;
IN	count	-	максимальное число принимаемых элементов;
IN	datatype	-	тип элементов принимаемого сообщения;
IN	source	-	номер процесса-отправителя;
IN	tag	-	идентификатор сообщения;
IN	comm	-	коммуникатор области связи;
OUT	status	-	атрибуты принятого сообщения.


int MPI_Barrier(MPI_Comm comm)
Функция барьерной синхронизации MPI_BARRIER блокирует вызывающий процесс, пока все процессы группы не вызовут её.
В каждом процессе управление возвращается только тогда, когда все процессы в группе вызовут процедуру.

Широковещательная рассылка данных выполняется с помощью функции MPI_Bcast. Процесс с номером root рассылает сообщение
из своего буфера передачи всем процессам области связи коммуникатора comm.

Математические операции над блоками данных, распределенных по процессорам, называют глобальными операциями редукции.
В общем случае операцией редукции называется операция, аргументом которой является вектор,
а результатом - скалярная величина, полученная применением некоторой математической операции ко всем компонентам вектора

В MPI глобальные операции редукции представлены в нескольких вариантах:
с сохранением результата в адресном пространстве одного процесса (MPI_Reduce).
с сохранением результата в адресном пространстве всех процессов (MPI_Allreduce).