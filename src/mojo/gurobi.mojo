
from memory import ArcPointer, UnsafePointer
from collections.string import String
from sys.ffi import DLHandle, OpaquePointer, external_call

var handle = DLHandle('libgurobi.so')


struct Environment(Movable):
    var ptr: OpaquePointer

    fn __init__(out self) raises:
        self.ptr = OpaquePointer()
        var errno = handle.call[
            'GRBloadenv',
            Int,
            UnsafePointer[OpaquePointer],
            UnsafePointer[UInt8]
        ](
            UnsafePointer.address_of(self.ptr),
            UnsafePointer[UInt8]()
        )
        if errno != 0:
            raise Error('Error creating Gurobi environment')

    fn __moveinit__(out self, owned existing: Self):
        self.ptr = existing.ptr

    fn __del__(owned self):
        _ = handle.call['GRBfreeenv', Int, OpaquePointer](self.ptr)
    
    @staticmethod
    fn create() raises -> ArcPointer[Environment]:
        """Creates a new environment."""
        return ArcPointer(Environment())

    # parameter management

    fn _get_int_param[name: StringLiteral](self) raises -> Int:
        var result: Int32 = 0
        var errno = handle.call[
            'GRBgetintparam',
            Int,
            OpaquePointer,
            UnsafePointer[Int8],
            UnsafePointer[Int32],
        ](
            self.ptr,
            name.unsafe_cstr_ptr(),
            UnsafePointer[Int32].address_of(result),
        )
        if errno != 0:
            raise Error('Error getting int parameter \'' + name + '\': ' + String(errno))
        return Int(result)
    
    fn _set_int_param[name: StringLiteral](self, val: Int) raises:
        var errno = handle.call[
            'GRBsetintparam',
            Int,
            OpaquePointer,
            UnsafePointer[Int8],
            Int32,
        ](
            self.ptr,
            name.unsafe_cstr_ptr(),
            val,
        )
        if errno != 0:
            raise Error('Error setting int parameter \'' + name + '\': ' + String(errno))

    fn get_threads(self) raises -> Int:
        return self._get_int_param['Threads']()

    fn set_threads(self, val: Int) raises:
        self._set_int_param['Threads'](val)
    
    # method
    # -1: automatic
    #  0: primal simplex
    #  1: dual simplex
    #  2: barrier
    #  3: concurrent
    #  4: deterministic concurrent
    #  5: deterministic concurrent simplex (deprecated)

    fn get_method(self) raises -> Int:
        return self._get_int_param['Method']()

    fn set_method(self, val: Int) raises:
        self._set_int_param['Method'](val)
    
    # output flag
    #  0: output disable
    #  1: output enabled

    fn get_output_flag(self) raises -> Int:
        return self._get_int_param['OutputFlag']()

    fn set_output_flag(self, val: Int) raises:
        self._set_int_param['OutputFlag'](val)


struct Model(Movable):
    var env: ArcPointer[Environment]
    var ptr: OpaquePointer

    fn __init__(out self, env: ArcPointer[Environment], read name: String) raises:
        self.env = env
        self.ptr = OpaquePointer()
        var errno = handle.call[
            'GRBnewmodel',
            Int,
            OpaquePointer,
            UnsafePointer[OpaquePointer],
            UnsafePointer[Int8],
            Int,
            UnsafePointer[Float64],
            UnsafePointer[Float64],
            UnsafePointer[Float64],
            UnsafePointer[Int32],
            UnsafePointer[UnsafePointer[Int8]]
        ](
            env[].ptr,
            UnsafePointer.address_of(self.ptr),
            name.unsafe_cstr_ptr(),
            0,
            UnsafePointer[Float64](),
            UnsafePointer[Float64](),
            UnsafePointer[Float64](),
            UnsafePointer[Int32](),
            UnsafePointer[UnsafePointer[Int8]]()
        )
        if errno != 0:
            raise Error('Error creating Gurobi model: ' + String(errno))

    fn __moveinit__(out self, owned existing: Self):
        self.env = existing.env^
        self.ptr = existing.ptr

    fn __del__(owned self):
        _ = handle.call['GRBfreemodel', Int, OpaquePointer](self.ptr)
    
    @staticmethod
    fn create(env: ArcPointer[Environment], read name: String) raises -> ArcPointer[Model]:
        """Creates a new environment."""
        return ArcPointer(Model(env, name))

    fn add_var(self, nnz: Int, idx: List[Int32], val: List[Float64], obj: Float64, lb: Float64, ub: Float64, vtype: Int, name: String) raises:
        # mod.add_var(idx=idx, val=val, obj=1.0, lb=0.0, ub=Float64('inf'), vtype=ord('C'), name=f"x{j:04d}")
        var errno = handle.call[
            'GRBaddvar',
            Int,
            OpaquePointer,
            Int,
            UnsafePointer[Int32],
            UnsafePointer[Float64],
            Float64,
            Float64,
            Float64,
            Int8,
            UnsafePointer[Int8]
        ](
            self.ptr,
            nnz,
            idx.unsafe_ptr(),
            val.unsafe_ptr(),
            obj,
            lb,
            ub,
            vtype,
            name.unsafe_cstr_ptr()
        )
        if errno != 0:
            raise Error('Error adding variable: ' + String(errno))

    fn add_constr(self, sense: Int, rhs: Float64, read name: String) raises:
        var errno = handle.call[
            'GRBaddconstr',
            Int,
            OpaquePointer,
            Int,
            UnsafePointer[Int32],
            UnsafePointer[Float64],
            Int8,
            Float64,
            UnsafePointer[Int8],
        ](
            self.ptr,
            0,
            UnsafePointer[Int32](),
            UnsafePointer[Float64](),
            sense,
            rhs,
            name.unsafe_cstr_ptr()
        )
        if errno != 0:
            raise Error('Error adding constraint: ' + String(errno))

    fn update(self) raises:
        var errno = handle.call[
            'GRBupdatemodel',
            Int,
            OpaquePointer,
        ](
            self.ptr
        )
        if errno != 0:
            raise Error('Error updating model: ' + String(errno))

    fn optimize(self) raises:
        var errno = handle.call[
            'GRBoptimize',
            Int,
            OpaquePointer,
        ](
            self.ptr
        )
        if errno != 0:
            raise Error('Error optimizing model: ' + String(errno))

    fn write(self, read filename: String) raises:
        var errno = handle.call[
            'GRBwrite',
            Int,
            OpaquePointer,
            UnsafePointer[Int8],
        ](
            self.ptr,
            filename.unsafe_cstr_ptr(),
        )
        if errno != 0:
            raise Error('Error writing model: ' + String(errno))

    # attribute utilities

    fn _get_int_attr[name: StringLiteral](self) raises -> Int:
        var result: Int32 = 0
        var errno = handle.call[
            'GRBgetintattr',
            Int,
            OpaquePointer,
            UnsafePointer[Int8],
            UnsafePointer[Int32],
        ](
            self.ptr,
            name.unsafe_cstr_ptr(),
            UnsafePointer[Int32].address_of(result),
        )
        if errno != 0:
            raise Error('Error getting int attribute \'' + name + '\': ' + String(errno))
        return Int(result)
    
    fn _set_int_attr[name: StringLiteral](self, val: Int) raises:
        var errno = handle.call[
            'GRBsetintattr',
            Int,
            OpaquePointer,
            UnsafePointer[Int8],
            Int32,
        ](
            self.ptr,
            name.unsafe_cstr_ptr(),
            val,
        )
        if errno != 0:
            raise Error('Error setting int attribute \'' + name + '\': ' + String(errno))

    fn _get_dbl_attr[name: StringLiteral](self) raises -> Float64:
        var result: Float64 = 0.0
        var errno = handle.call[
            'GRBgetdblattr',
            Int,
            OpaquePointer,
            UnsafePointer[Int8],
            UnsafePointer[Float64],
        ](
            self.ptr,
            name.unsafe_cstr_ptr(),
            UnsafePointer[Float64].address_of(result),
        )
        if errno != 0:
            raise Error('Error getting double attribute \'' + name + '\': ' + String(errno))
        return result
    
    fn _set_dbl_attr[name: StringLiteral](self, val: Float64) raises:
        var errno = handle.call[
            'GRBsetdblattr',
            Int,
            OpaquePointer,
            UnsafePointer[Int8],
            Float64,
        ](
            self.ptr,
            name.unsafe_cstr_ptr(),
            val,
        )
        if errno != 0:
            raise Error('Error setting double attribute \'' + name + '\': ' + String(errno))

    fn _get_int_attr_element[name: StringLiteral](self, i: Int) raises -> Int:
        var result: Int32 = 0
        var errno = handle.call[
            'GRBgetintattrelement',
            Int,
            OpaquePointer,
            UnsafePointer[Int8],
            Int32,
            UnsafePointer[Int32],
        ](
            self.ptr,
            name.unsafe_cstr_ptr(),
            i,
            UnsafePointer[Int32].address_of(result),
        )
        if errno != 0:
            raise Error('Error getting int attribute \'' + name + '\': ' + String(errno))
        return Int(result)
    
    fn _set_int_attr_element[name: StringLiteral](self, i: Int, val: Int) raises:
        var errno = handle.call[
            'GRBsetintattrelement',
            Int,
            OpaquePointer,
            UnsafePointer[Int8],
            Int32,
            Int32,
        ](
            self.ptr,
            name.unsafe_cstr_ptr(),
            i,
            val,
        )
        if errno != 0:
            raise Error('Error setting int attribute \'' + name + '\': ' + String(errno))

    fn _get_dbl_attr_element[name: StringLiteral](self, i: Int) raises -> Float64:
        var result: Float64 = 0.0
        var errno = handle.call[
            'GRBgetdblattrelement',
            Int,
            OpaquePointer,
            UnsafePointer[Int8],
            Int32,
            UnsafePointer[Float64],
        ](
            self.ptr,
            name.unsafe_cstr_ptr(),
            i,
            UnsafePointer[Float64].address_of(result),
        )
        if errno != 0:
            raise Error('Error getting double attribute \'' + name + '\': ' + String(errno))
        return result

    # model attributes

    # Status
    #  1: LOADED
    #  2: OPTIMAL
    #  3: INFEASIBLE
    #  4: INF_OR_UNBD
    #  5: UNBOUNDED
    #  6: CUTOFF
    #  7: ITERATION_LIMIT
    #  8: NODE_LIMIT
    #  9: TIME_LIMIT
    # 10: SOLUTION_LIMIT
    # 11: INTERRUPTED
    # 12: NUMERIC
    # 13: SUBOPTIMAL
    # 14: INPROGRESS
    # 15: USER_OBJ_LIMIT
    # 16: WORK_LIMIT
    # 17: MEM_LIMIT
    fn get_status(self) raises -> Int:
        return self._get_int_attr['Status']()

    fn get_model_sense(self) raises -> Int:
        return self._get_int_attr['ModelSense']()

    fn set_model_sense(self, sense: Int) raises:
        self._set_int_attr['ModelSense'](sense)

    fn get_objval(self) raises -> Float64:
        return self._get_dbl_attr['ObjVal']()

    # variable attributes

    fn get_lb(self, i: Int) raises -> Float64:
        return self._get_dbl_attr_element['LB'](i)

    fn get_ub(self, i: Int) raises -> Float64:
        return self._get_dbl_attr_element['UB'](i)

    fn get_obj(self, i: Int) raises -> Float64:
        return self._get_dbl_attr_element['Obj'](i)

    fn get_x(self, i: Int) raises -> Float64:
        return self._get_dbl_attr_element['X'](i)

    # VBasis
    #  0: basic
    # -1: non-basic

    fn get_vbasis(self, i: Int) raises -> Int:
        return self._get_int_attr_element['VBasis'](i)

    # constraint attributes

    fn get_sense(self, i: Int) raises -> Int8:
        return self._get_int_attr_element['Sense'](i)

    fn get_rhs(self, i: Int) raises -> Float64:
        return self._get_dbl_attr_element['RHS'](i)

    fn get_pi(self, i: Int) raises -> Float64:
        return self._get_dbl_attr_element['Pi'](i)

    # CBasis
    #  0: basic
    # -1: non-basic at lower bound
    # -2: non-basic at upper bound
    # -3: super-basic

    fn get_cbasis(self, i: Int) raises -> Int:
        return self._get_int_attr_element['CBasis'](i)



