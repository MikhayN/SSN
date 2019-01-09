module PQ 

/// Create binary tree with MAX on top
type MaxIndexPriorityQueue<'T when 'T : comparison>(n:int) = 

    let objects : 'T []      = Array.zeroCreate n //ResizeArray<'T>(n)
    let heap    : int []     = Array.zeroCreate (n+1) //ResizeArray<int>(n)
    let heapInverse : int [] = Array.zeroCreate n//ResizeArray<int>(n)
    let mutable m_count      = 0

    let swapInplace i j =
        // swap elements in heap
        let tmp = heap.[i]
        heap.[i] <- heap.[j]
        heap.[j] <- tmp

        // reset inverses
        heapInverse.[heap.[i]] <- i
        heapInverse.[heap.[j]] <- j

    let parent heapIndex = heapIndex / 2
    let firstChild heapIndex = heapIndex * 2 
    let secondChild heapIndex = heapIndex * 2 + 1

    let sortHeapDownward heapIndex = // we are looking if the object is smaller then its children 
        let rec loop heapIndex =
            let childIndex = firstChild heapIndex
            if (childIndex <= m_count) then
                let child = // choose the biggest of children, if there both availbale 
                    if (childIndex < m_count) && ( objects.[heap.[childIndex + 1]] > objects.[heap.[childIndex]] ) then // change here
                        childIndex+1 
                    else 
                        childIndex
                // swap with child if the child is bigger
                if ( objects.[heap.[child]] > objects.[heap.[heapIndex]] ) then // change here
                    swapInplace child heapIndex
                    loop child
        loop heapIndex

    let sortHeapUpward heapIndex = // check if the object bigger than its parent, if yes, move it upper
        let rec loop heapIndex =
            let parentIndex = parent( heapIndex )
            if (heapIndex > 1 && 
                    objects.[heap.[heapIndex]] > objects.[heap.[parentIndex]] ) then // change here 
                // swap this node with its parent
                swapInplace heapIndex parentIndex
                // reset iterator to be at parents old position
                // (child's new position)
                loop parentIndex
        loop heapIndex

    let sortUpward realIndex   = sortHeapUpward   heapInverse.[realIndex]
    let sortDownward realIndex = sortHeapDownward heapInverse.[realIndex]
    
    /// Clear the counter, get ready to refill the queue
    member this.Clear() = m_count <- 0

    /// Increase the value at the current index
    member this.IncreaseValueAtIndex index (obj : 'T) =
        if not (index < objects.Length && index >= 0) then 
            failwithf "IndexedPriorityQueue.DecreaseIndex: Index %i out of range" index
        if (obj <= objects.[index]) then // check if new value really bigger
            failwithf "IndexedPriorityQueue.DecreaseIndex: object '%A' isn't greater than current value '%A'" obj objects.[index]
        objects.[index] <- obj
        sortUpward index

    /// Decrease the value at the current index
    member this.DecreaseValueAtIndex index (obj:'T) =
        if not (index < objects.Length && index >= 0) then 
            failwithf "IndexedPriorityQueue.DecreaseIndex: Index %i out of range" index
        if (obj >= objects.[index]) then // check if new value really smaller
            failwithf "IndexedPriorityQueue.DecreaseIndex: object '%A' isn't less than current value '%A'" obj objects.[index]
        objects.[index] <- obj
        sortDownward index

    /// Updates the value at the given index. Note that this function is not
    /// as efficient as the DecreaseIndex/IncreaseIndex methods, but is
    /// best when the value at the index is not known
    member this.Set index (obj:'T) =
        if ( obj >= objects.[index] ) then // change here if change MAX MIN sorting
            this.IncreaseValueAtIndex index obj 
        else 
            this.DecreaseValueAtIndex index obj

    /// Removes the top element from the queue 
    member this.Pop () =
        //if (m_count > 0) then failwith "IndexedPriorityQueue.Pop: Queue is empty"
        if (m_count = 0) then 
            Unchecked.defaultof<'T>
        else
            // swap front to back for removal
            swapInplace 1 (m_count)
            m_count <- m_count - 1
            // re-sort heap
            sortHeapDownward 1
            // return popped object
            //m_objects[m_heap[m_count + 1]]
            objects.[heap.[m_count+1]]

    /// Removes the element with given index from the queue 
    member this.Remove index =
        //if (m_count > 0) then failwith "IndexedPriorityQueue.Pop: Queue is empty"
        if ((heapInverse.[index]) > m_count) then 
            failwith "IndexedPriorityQueue.Remove: The element was already removed"
        else
            // swap front to back for removal
            swapInplace heapInverse.[index] (m_count)
            m_count <- m_count - 1
            // re-sort heap
            sortHeapDownward heapInverse.[index]

    /// Removes the element with given index from the queue 
    member this.TryRemove index =
        //if (m_count > 0) then failwith "IndexedPriorityQueue.Pop: Queue is empty"
        if ((heapInverse.[index]) <= m_count) then 
            
            // swap front to back for removal
            swapInplace heapInverse.[index] (m_count)
            m_count <- m_count - 1
            // re-sort heap
            sortHeapDownward heapInverse.[index]
    
    /// Put indexOut element outside of the heap and move indexIn element inside the heap
    member this.Swap indexOut indexIn =
        if ((heapInverse.[indexOut]) <= m_count) then 
            failwith "IndexedPriorityQueue.Return: The element was already removed from heap"
        elif ((heapInverse.[indexIn]) > m_count) then 
            failwith "IndexedPriorityQueue.Return: The element was already inside heap"  
        else   
            swapInplace (heapInverse.[indexOut]) (heapInverse.[indexIn])
            if ( objects.[indexIn] > objects.[indexOut] ) then // change here if change MAX MIN sorting
                sortUpward indexOut
            else 
                sortDownward indexOut

    /// Removes the element with given index from the queue 
    member this.RemoveGroup (indeces: int []) =
        let heapIDs = 
            [|for id in indeces ->
                let c = heapInverse.[id]
                if (c > m_count) then 
                    failwith "IndexedPriorityQueue.Remove: The element was already removed"
                else
                    // swap front to back for removal
                    swapInplace c (m_count)
                    m_count <- m_count - 1
                    c  |] 
            |> Array.sort
        // re-sort heap
        for id in heapIDs do
            sortHeapUpward id

    /// Removes the element with given index from the queue 
    member this.TryRemoveGroup (indeces: int []) =
        //let heapIDs = 
        //    [|for id in indeces ->
        //        let c = heapInverse.[id]
        //        if (c > m_count) then 
        //            failwith "IndexedPriorityQueue.Remove: The element was already removed"
        //        else
        //            // swap front to back for removal
        //            swapInplace c (m_count)
        //            m_count <- m_count - 1
        //            c  |] 
        //    |> Array.sort
        let heapIDs =
            indeces
            |> Array.map (fun id -> heapInverse.[id])
            |> Array.filter (fun c -> c <= m_count)
            |> Array.map (fun c -> 
                swapInplace c (m_count)
                m_count <- m_count - 1
                c)
            |> Array.sort
        // re-sort heap
        for id in heapIDs do
            sortHeapUpward id

    member this.ReturnGroup (indeces: int[]) =
        let heapIDs = 
            [|for id in indeces ->
                let c = heapInverse.[id]
                if (c <= m_count) then 
                    failwith "IndexedPriorityQueue.ReturnGroup: The element was already inside"
                else
                    // swap front to back for removal
                    swapInplace c (m_count+1)
                    m_count <- m_count + 1
                    c  |] 
            |> Array.sort
        // re-sort heap
        for id in heapIDs do
            sortHeapUpward id

    member this.TryReturnGroup (indeces: int[]) =
        let heapIDs =
            indeces
            |> Array.map (fun id -> heapInverse.[id])
            |> Array.filter (fun c -> c > m_count)
            |> Array.map (fun c -> 
                swapInplace c (m_count+1)
                m_count <- m_count + 1
                c)
            |> Array.sort
        // re-sort heap
        for id in heapIDs do
            sortHeapUpward id

    /// Removes all elements except those with given index from the queue 
    member this.LeaveGroup (indeces: int []) =
        let x = indeces.Length
        indeces
        |> Array.map (fun i -> heapInverse.[i])
        |> Array.sort
        |> Array.iteri (fun i c ->
            //if (c > m_count) then 
            //    failwith "IndexedPriorityQueue.Remove: The element was already removed"
            //else
                // swap front to back for removal
                swapInplace c (i+1) )
        // re-sort heap
        m_count <- x
        for id=2 to x do
            sortHeapUpward id

    /// Gets the top element of the queue
    member this.Top() = 
        // top of heap [first element is 1, not 0]
        objects.[heap.[1]]

    /// Gets the index of the top element of the queue
    member this.TopIndex() = 
        // top of heap [first element is 1, not 0]
        heap.[1]

    /// Inserts a new value with the given index in the queue
    member this.Insert index (value:'T) =
        //if (index < objects.Length && index >= 0) then failwithf "IndexedPriorityQueue.Insert: Index %i out of range" index
        m_count <- m_count + 1
        
        // add object
        objects.[index] <- value

        // add to heap
        heapInverse.[index] <- m_count
        heap.[m_count] <- index

        // update heap
        sortHeapUpward m_count // the new item is by default in the bottom of the tree, as the smallest. So check, if it is not the case, then move it up
        
    member this.Length = 
        m_count

    /// Returns an item of the set with given index
    member this.Item 
        with get index = 
            if not (index < objects.Length && index >= 0) then failwithf "IndexedPriorityQueue.Item: Index %i out of range" index
            objects.[index]
        and  set index value = 
            if not (index < objects.Length && index >= 0) then failwithf "IndexedPriorityQueue.Item: Index %i out of range" index
            this.Set index  value 
    
    /// Returns an item of the heap with given index
    member this.HeapItem index =
        if (index <= 0) || (index > m_count) then failwithf "IndexedPriorityQueue.Item: Index %i out of range" index
        objects.[heap.[index]]

    /// Returns an index in original objects order of the heap element with given index
    member this.HeapItemIndex index =
        if (index <= 0) || (index > m_count) then failwithf "IndexedPriorityQueue.Item: Index %i out of range" index
        heap.[index]

    /// Set internal properties of the type directly    
    member internal this.SetData (_objects: 'T []) (_heap:int[]) (_heapInverse:int[]) (_m_count: int) = 
        for i=0 to (n-1) do
            objects.[i]       <- _objects.[i] 
            heap.[i]          <- _heap.[i]
            heapInverse.[i]   <- _heapInverse.[i]
        heap.[objects.Length] <- _heap.[objects.Length]
        m_count <- _m_count
    
    /// deep Copy the variable into a new one
    member this.DeepCopy() = 
        let temp = new MaxIndexPriorityQueue<'T>(n)
        temp.SetData objects heap heapInverse m_count
        temp

/// Create binary tree with MAX on top
type MinIndexPriorityQueue<'T when 'T : comparison>(n:int) = 

    let objects : 'T []      = Array.zeroCreate n //ResizeArray<'T>(n)
    let heap    : int []     = Array.zeroCreate (n+1) //ResizeArray<int>(n)
    let heapInverse : int [] = Array.zeroCreate n//ResizeArray<int>(n)
    let mutable m_count      = 0

    let swapInplace i j =
        // swap elements in heap
        let tmp = heap.[i]
        heap.[i] <- heap.[j]
        heap.[j] <- tmp

        // reset inverses
        heapInverse.[heap.[i]] <- i
        heapInverse.[heap.[j]] <- j

    let parent heapIndex = heapIndex / 2
    let firstChild heapIndex = heapIndex * 2 
    let secondChild heapIndex = heapIndex * 2 + 1

    let sortHeapDownward heapIndex = // we are looking if the object is bigger then its children 
        let rec loop heapIndex =
            let childIndex = firstChild heapIndex
            if (childIndex <= m_count) then
                let child = // choose the smallest of children, if there both availbale 
                    if (childIndex < m_count) && ( objects.[heap.[childIndex + 1]] < objects.[heap.[childIndex]] ) then // change here
                        childIndex+1 
                    else 
                        childIndex
                // swap with child if the child is smaller
                if ( objects.[heap.[child]] < objects.[heap.[heapIndex]] ) then // change here
                    swapInplace child heapIndex
                    loop child
        loop heapIndex

    let sortHeapUpward heapIndex = // check if the object smaller than its parent, if yes, move it upper
        let rec loop heapIndex =
            let parentIndex = parent( heapIndex )
            if (heapIndex > 1 && 
                    objects.[heap.[heapIndex]] < objects.[heap.[parentIndex]] ) then // change here 
                // swap this node with its parent
                swapInplace heapIndex parentIndex
                // reset iterator to be at parents old position
                // (child's new position)
                loop parentIndex
        loop heapIndex

    let sortUpward realIndex   = sortHeapUpward   heapInverse.[realIndex]
    let sortDownward realIndex = sortHeapDownward heapInverse.[realIndex]
    
    /// Clear the counter, get ready to refill the queue
    member this.Clear() = m_count <- 0
    
    /// Increase the value at the current index
    member this.IncreaseValueAtIndex index (obj : 'T) =
        if not (index < objects.Length && index >= 0) then 
            failwithf "IndexedPriorityQueue.DecreaseIndex: Index %i out of range" index
        if (obj <= objects.[index]) then // check if new value really bigger
            failwithf "IndexedPriorityQueue.DecreaseIndex: object '%A' isn't greater than current value '%A'" obj objects.[index]
        objects.[index] <- obj
        sortDownward index

    /// Decrease the value at the current index
    member this.DecreaseValueAtIndex index (obj:'T) =
        if not (index < objects.Length && index >= 0) then 
            failwithf "IndexedPriorityQueue.DecreaseIndex: Index %i out of range" index
        if (obj >= objects.[index]) then // check if new value really smaller
            failwithf "IndexedPriorityQueue.DecreaseIndex: object '%A' isn't less than current value '%A'" obj objects.[index]
        objects.[index] <- obj
        sortUpward index

    /// Updates the value at the given index. Note that this function is not
    /// as efficient as the DecreaseIndex/IncreaseIndex methods, but is
    /// best when the value at the index is not known
    member this.Set index (obj:'T) =
        if ( obj >= objects.[index] ) then // change here if change MAX MIN sorting
            this.IncreaseValueAtIndex index obj 
        else 
            this.DecreaseValueAtIndex index obj

    /// Removes the top element from the queue 
    member this.Pop () =
        //if (m_count > 0) then failwith "IndexedPriorityQueue.Pop: Queue is empty"
        if (m_count = 0) then 
            Unchecked.defaultof<'T>
        else
            // swap front to back for removal
            swapInplace 1 (m_count)
            m_count <- m_count - 1
            // re-sort heap
            sortHeapDownward 1
            // return popped object
            //m_objects[m_heap[m_count + 1]]
            objects.[heap.[m_count+1]]

    /// Removes the element with given index from the queue 
    member this.Remove index =
        //if (m_count > 0) then failwith "IndexedPriorityQueue.Pop: Queue is empty"
        if ((heapInverse.[index]) > m_count) then 
            failwith "IndexedPriorityQueue.Remove: The element was already removed"
        else
            // swap front to back for removal
            swapInplace heapInverse.[index] (m_count)
            m_count <- m_count - 1
            // re-sort heap
            sortHeapDownward heapInverse.[index]

    /// Removes the element with given index from the queue 
    member this.TryRemove index =
        //if (m_count > 0) then failwith "IndexedPriorityQueue.Pop: Queue is empty"
        if ((heapInverse.[index]) <= m_count) then 
            
            // swap front to back for removal
            swapInplace heapInverse.[index] (m_count)
            m_count <- m_count - 1
            // re-sort heap
            sortHeapDownward heapInverse.[index]
    
    /// Put indexOut element outside of the heap and move indexIn element inside the heap
    member this.Swap indexOut indexIn =
        if ((heapInverse.[indexOut]) <= m_count) then 
            failwith "IndexedPriorityQueue.Return: The element was already removed from heap"
        elif ((heapInverse.[indexIn]) > m_count) then 
            failwith "IndexedPriorityQueue.Return: The element was already inside heap"  
        else   
            swapInplace (heapInverse.[indexOut]) (heapInverse.[indexIn])
            if ( objects.[indexIn] > objects.[indexOut] ) then // change here if change MAX MIN sorting
                sortUpward indexOut
            else 
                sortDownward indexOut

    /// Removes the element with given index from the queue 
    member this.TryRemoveGroup (indeces: int []) =
        //let heapIDs = 
        //    [|for id in indeces ->
        //        let c = heapInverse.[id]
        //        if (c > m_count) then 
        //            failwith "IndexedPriorityQueue.Remove: The element was already removed"
        //        else
        //            // swap front to back for removal
        //            swapInplace c (m_count)
        //            m_count <- m_count - 1
        //            c  |] 
        //    |> Array.sort
        let heapIDs =
            indeces
            |> Array.map (fun id -> heapInverse.[id])
            |> Array.filter (fun c -> c <= m_count)
            |> Array.map (fun c -> 
                swapInplace c (m_count)
                m_count <- m_count - 1
                c)
            |> Array.sort
        // re-sort heap
        for id in heapIDs do
            sortHeapUpward id

    member this.ReturnGroup (indeces: int[]) =
        let heapIDs = 
            [|for id in indeces ->
                let c = heapInverse.[id]
                if (c <= m_count) then 
                    failwith "IndexedPriorityQueue.ReturnGroup: The element was already inside"
                else
                    // swap front to back for removal
                    swapInplace c (m_count+1)
                    m_count <- m_count + 1
                    c  |] 
            |> Array.sort
        // re-sort heap
        for id in heapIDs do
            sortHeapUpward id

    member this.TryReturnGroup (indeces: int[]) =
        let heapIDs =
            indeces
            |> Array.map (fun id -> heapInverse.[id])
            |> Array.filter (fun c -> c > m_count)
            |> Array.map (fun c -> 
                swapInplace c (m_count+1)
                m_count <- m_count + 1
                c)
            |> Array.sort
        // re-sort heap
        for id in heapIDs do
            sortHeapUpward id

    /// Removes all elements except those with given index from the queue 
    member this.LeaveGroup (indeces: int []) =
        let x = indeces.Length
        indeces
        |> Array.map (fun i -> heapInverse.[i])
        |> Array.sort
        |> Array.iteri (fun i c ->
            //if (c > m_count) then 
            //    failwith "IndexedPriorityQueue.Remove: The element was already removed"
            //else
                // swap front to back for removal
                swapInplace c (i+1) )
        // re-sort heap
        m_count <- x
        for id=2 to x do
            sortHeapUpward id

    /// Gets the top element of the queue
    member this.Top() = 
        // top of heap [first element is 1, not 0]
        objects.[heap.[1]]

    /// Gets the index of the top element of the queue
    member this.TopIndex() = 
        // top of heap [first element is 1, not 0]
        heap.[1]

    /// Inserts a new value with the given index in the queue
    member this.Insert index (value:'T) =
        //if (index < objects.Length && index >= 0) then failwithf "IndexedPriorityQueue.Insert: Index %i out of range" index
        m_count <- m_count + 1
        
        // add object
        objects.[index] <- value

        // add to heap
        heapInverse.[index] <- m_count
        heap.[m_count] <- index

        // update heap
        sortHeapUpward m_count // the new item is by default in the bottom of the tree, the biggest. So check, if it is not the case, then move it up
        
    member this.Length = 
        m_count

    /// Returns an item of the set with given index
    member this.Item 
        with get index = 
            if not (index < objects.Length && index >= 0) then failwithf "IndexedPriorityQueue.Item: Index %i out of range" index
            objects.[index]
        and  set index value = 
            if not (index < objects.Length && index >= 0) then failwithf "IndexedPriorityQueue.Item: Index %i out of range" index
            this.Set index  value 
    
    /// Returns an item of the heap with given index
    member this.HeapItem index =
        if (index <= 0) || (index > m_count) then failwithf "IndexedPriorityQueue.Item: Index %i out of range" index
        objects.[heap.[index]]

    /// Returns an index in original objects order of the heap element with given index
    member this.HeapItemIndex index =
        if (index <= 0) || (index > m_count) then failwithf "IndexedPriorityQueue.Item: Index %i out of range" index
        heap.[index]

    /// Set internal properties of the type directly    
    member internal this.SetData (_objects: 'T []) (_heap:int[]) (_heapInverse:int[]) (_m_count: int) = 
        for i=0 to (n-1) do
            objects.[i]       <- _objects.[i] 
            heap.[i]          <- _heap.[i]
            heapInverse.[i]   <- _heapInverse.[i]
        heap.[objects.Length] <- _heap.[objects.Length]
        m_count <- _m_count
    
    /// deep Copy the variable into a new one
    member this.DeepCopy() = 
        let temp = new MinIndexPriorityQueue<'T>(n)
        temp.SetData objects heap heapInverse m_count
        temp