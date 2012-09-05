<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="text" omit-xml-declaration="yes"/>

<xsl:variable name="std_namespace_id" select="/GCC_XML/Namespace[@name='std']/@id"/>
<xsl:variable name="openmm_namespace_id" select="/GCC_XML/Namespace[@name='OpenMM']/@id"/>
<xsl:variable name="void_type_id" select="/GCC_XML/FundamentalType[@name='void']/@id"/>
<xsl:variable name="int_type_id" select="/GCC_XML/FundamentalType[@name='int']/@id"/>
<xsl:variable name="double_type_id" select="/GCC_XML/FundamentalType[@name='double']/@id"/>
<xsl:variable name="bool_type_id" select="/GCC_XML/FundamentalType[@name='bool']/@id"/>
<xsl:variable name="char_type_id" select="/GCC_XML/FundamentalType[@name='char']/@id"/>
<xsl:variable name="string_type_id" select="/GCC_XML/*[@name='string' and @context=$std_namespace_id]/@id"/>
<xsl:variable name="vec3_type_id" select="/GCC_XML/*[@name='Vec3' and @context=$openmm_namespace_id]/@id"/>
<xsl:variable name="const_char_type_id" select="/GCC_XML/CvQualifiedType[@type=$char_type_id]/@id"/>
<xsl:variable name="ptr_const_char_type_id" select="/GCC_XML/PointerType[@type=$const_char_type_id]/@id"/>
<xsl:variable name="vector_string_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;std::basic_string')]/@id"/>
<xsl:variable name="vector_vec3_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;OpenMM::Vec3')]/@id"/>
<xsl:variable name="vector_bond_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;std::pair&lt;int, int')]/@id"/>
<xsl:variable name="map_parameter_type_id" select="/GCC_XML/Class[starts-with(@name, 'map&lt;std::basic_string') and contains(@name, 'double')]/@id"/>
<xsl:variable name="map_property_type_id" select="/GCC_XML/Class[starts-with(@name, 'map&lt;std::basic_string') and not(contains(@name, 'double'))]/@id"/>
<xsl:variable name="vector_double_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;double')]/@id"/>
<xsl:variable name="vector_int_type_id" select="/GCC_XML/Class[starts-with(@name, 'vector&lt;int')]/@id"/>
<xsl:variable name="newline">
<xsl:text>
</xsl:text>
</xsl:variable>

<!-- Do not generate functions for the following classes -->
<xsl:variable name="skip_classes" select="('Vec3', 'Kernel', 'Stream', 'KernelImpl', 'StreamImpl', 'KernelFactory', 'StreamFactory')"/>
<!-- Do not generate the following functions -->
<xsl:variable name="skip_methods" select="('OpenMM_Context_getState', 'OpenMM_Platform_loadPluginsFromDirectory', 'OpenMM_Context_createCheckpoint', 'OpenMM_Context_loadCheckpoint')"/>
<!-- Suppress any function which references any of the following classes -->
<xsl:variable name="hide_classes" select="('Kernel', 'Stream', 'KernelImpl', 'StreamImpl', 'KernelFactory', 'StreamFactory', 'ContextImpl')"/>

<!-- Main loop over all classes in the OpenMM namespace -->
<xsl:template match="/GCC_XML">
MODULE OpenMM_Types
    implicit none

    ! Global Constants
 <xsl:for-each select="Variable[@context=$openmm_namespace_id]">
    real*8, parameter :: OpenMM_<xsl:value-of select="concat(@name, ' = ', number(@init))"/>
 </xsl:for-each>

    ! Type Declarations
 <xsl:for-each select="(Class | Struct)[@context=$openmm_namespace_id and empty(index-of($skip_classes, @name))]">
    type OpenMM_<xsl:value-of select="@name"/>
        integer*8 :: handle = 0
    end type
 </xsl:for-each>
    type OpenMM_Vec3Array
        integer*8 :: handle = 0
    end type

    type OpenMM_StringArray
        integer*8 :: handle = 0
    end type

    type OpenMM_BondArray
        integer*8 :: handle = 0
    end type

    type OpenMM_ParameterArray
        integer*8 :: handle = 0
    end type

    type OpenMM_PropertyArray
        integer*8 :: handle = 0
    end type

    type OpenMM_DoubleArray
        integer*8 :: handle = 0
    end type

    type OpenMM_IntArray
        integer*8 :: handle = 0
    end type

    ! Enumerations

    integer*4, parameter :: OpenMM_False = 0
    integer*4, parameter :: OpenMM_True = 1
 <xsl:for-each select="Class[@context=$openmm_namespace_id and empty(index-of($skip_classes, @name))]">
  <xsl:variable name="class_id" select="@id"/>
  <xsl:variable name="class_name" select="@name"/>
  <xsl:for-each select="/GCC_XML/Enumeration[@context=$class_id and @access='public']">
   <xsl:call-template name="enumeration">
    <xsl:with-param name="class_name" select="$class_name"/>
   </xsl:call-template>
  </xsl:for-each>
 </xsl:for-each>
END MODULE OpenMM_Types

MODULE OpenMM
    use OpenMM_Types; implicit none
    interface

        ! OpenMM_Vec3
        subroutine OpenMM_Vec3_scale(vec, scale, result)
            use OpenMM_Types; implicit none
            real*8 vec(3)
            real*8 scale
            real*8 result(3)
        end subroutine

        ! OpenMM_Vec3Array
        subroutine OpenMM_Vec3Array_create(result, size)
            use OpenMM_Types; implicit none
            integer*4 size
            type (OpenMM_Vec3Array) result
        end subroutine
        subroutine OpenMM_Vec3Array_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) destroy
        end subroutine
        function OpenMM_Vec3Array_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            integer*4 OpenMM_Vec3Array_getSize
        end function
        subroutine OpenMM_Vec3Array_resize(target, size)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            integer*4 size
        end subroutine
        subroutine OpenMM_Vec3Array_append(target, vec)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            real*8 vec(3)
        end subroutine
        subroutine OpenMM_Vec3Array_set(target, index, vec)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            integer*4 index
            real*8 vec(3)
        end subroutine
        subroutine OpenMM_Vec3Array_get(target, index, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            integer*4 index
            real*8 result(3)
        end subroutine

        ! OpenMM_StringArray
        subroutine OpenMM_StringArray_create(result, size)
            use OpenMM_Types; implicit none
            integer*4 size
            type (OpenMM_StringArray) result
        end subroutine
        subroutine OpenMM_StringArray_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) destroy
        end subroutine
        function OpenMM_StringArray_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            integer*4 OpenMM_StringArray_getSize
        end function
        subroutine OpenMM_StringArray_resize(target, size)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            integer*4 size
        end subroutine
        subroutine OpenMM_StringArray_append(target, str)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            character(*) str
        end subroutine
        subroutine OpenMM_StringArray_set(target, index, str)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            integer*4 index
            character(*) str
        end subroutine
        subroutine OpenMM_StringArray_get(target, index, result)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            integer*4 index
            character(*) result
        end subroutine

        ! OpenMM_BondArray
        subroutine OpenMM_BondArray_create(result, size)
            use OpenMM_Types; implicit none
            integer*4 size
            type (OpenMM_BondArray) result
        end subroutine
        subroutine OpenMM_BondArray_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) destroy
        end subroutine
        function OpenMM_BondArray_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 OpenMM_BondArray_getSize
        end function
        subroutine OpenMM_BondArray_resize(target, size)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 size
        end subroutine
        subroutine OpenMM_BondArray_append(target, particle1, particle2)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 particle1
            integer*4 particle2
        end subroutine
        subroutine OpenMM_BondArray_set(target, index, particle1, particle2)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
        end subroutine
        subroutine OpenMM_BondArray_get(target, index, particle1, particle2)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
        end subroutine

        ! OpenMM_ParameterArray
        function OpenMM_ParameterArray_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_ParameterArray) target
            integer*4 OpenMM_ParameterArray_getSize
        end function
        subroutine OpenMM_ParameterArray_get(target, name, result)
            use OpenMM_Types; implicit none
            type (OpenMM_ParameterArray) target
            character(*) name
            character(*) result
        end subroutine

        ! OpenMM_PropertyArray
        function OpenMM_PropertyArray_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_ParameterArray) target
            integer*4 OpenMM_PropertyArray_getSize
        end function
        subroutine OpenMM_PropertyArray_get(target, name, result)
            use OpenMM_Types; implicit none
            type (OpenMM_PropertyArray) target
            character(*) name
            character(*) result
        end subroutine

<xsl:call-template name="primitive_array">
 <xsl:with-param name="element_type" select="'real*8'"/>
 <xsl:with-param name="name" select="'OpenMM_DoubleArray'"/>
</xsl:call-template>

<xsl:call-template name="primitive_array">
 <xsl:with-param name="element_type" select="'integer*4'"/>
 <xsl:with-param name="name" select="'OpenMM_IntArray'"/>
</xsl:call-template>

        ! These methods need to be handled specially, since their C++ APIs cannot be directly translated to Fortran.
        ! Unlike the C++ versions, the return value is allocated on the heap, and you must delete it yourself.
        subroutine OpenMM_Context_getState(target, types, enforcePeriodicBox, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
            integer*4 types
            integer*4 enforcePeriodicBox
            type (OpenMM_State) result
        end subroutine

        subroutine OpenMM_Platform_loadPluginsFromDirectory(directory, result)
            use OpenMM_Types; implicit none
            character(*) directory
            type (OpenMM_StringArray) result
        end subroutine

 <!-- Class members -->
 <xsl:for-each select="Class[@context=$openmm_namespace_id and empty(index-of($skip_classes, @name))]">
  <xsl:call-template name="class"/>
 </xsl:for-each>
    end interface
END MODULE OpenMM
</xsl:template>

<!-- Print out the declarations for a (Primitive)Array type -->
<xsl:template name="primitive_array">
 <xsl:param name="element_type"/>
 <xsl:param name="name"/>
        ! <xsl:value-of select="$name"/>
        subroutine <xsl:value-of select="$name"/>_create(result, size)
            use OpenMM_Types; implicit none
            integer*4 size
            type (<xsl:value-of select="$name"/>) result
        end subroutine
        subroutine <xsl:value-of select="$name"/>_destroy(destroy)
            use OpenMM_Types; implicit none
            type (<xsl:value-of select="$name"/>) destroy
        end subroutine
        function <xsl:value-of select="$name"/>_getSize(target)
            use OpenMM_Types; implicit none
            type (<xsl:value-of select="$name"/>) target
            integer*4 <xsl:value-of select="$name"/>_getSize
        end function
        subroutine <xsl:value-of select="$name"/>_resize(target, size)
            use OpenMM_Types; implicit none
            type (<xsl:value-of select="$name"/>) target
            integer*4 size
        end subroutine
        subroutine <xsl:value-of select="$name"/>_append(target, value)
            use OpenMM_Types; implicit none
            type (<xsl:value-of select="$name"/>) target
            <xsl:value-of select="$element_type"/> value
        end subroutine
        subroutine <xsl:value-of select="$name"/>_set(target, index, value)
            use OpenMM_Types; implicit none
            type (<xsl:value-of select="$name"/>) target
            integer*4 index
            <xsl:value-of select="$element_type"/> value
        end subroutine
        subroutine <xsl:value-of select="$name"/>_get(target, index, result)
            use OpenMM_Types; implicit none
            type (<xsl:value-of select="$name"/>) target
            integer*4 index
            <xsl:value-of select="$element_type"/> result
        end subroutine
</xsl:template>

<!-- Print out information for a class -->
<xsl:template name="class">
 <xsl:variable name="class_name" select="@name"/>
 <xsl:variable name="class_id" select="@id"/>

        ! OpenMM::<xsl:value-of select="@name"/>
 <!-- Constructors and destructor -->
 <xsl:if test="not(@abstract=1)">
  <xsl:variable name="constructors" select="/GCC_XML/Constructor[@context=$class_id and @access='public' and not(@artificial='1')]"/>
  <xsl:for-each select="$constructors">
   <xsl:call-template name="constructor">
    <xsl:with-param name="suffix" select="if (position() > 1) then concat('_', position()) else ''"/>
   </xsl:call-template>
  </xsl:for-each>
  <xsl:call-template name="destructor"/>
 </xsl:if>
 <!-- Methods -->
 <xsl:variable name="methods" select="/GCC_XML/Method[@context=$class_id and @access='public']"/>
 <xsl:for-each select="$methods">
  <xsl:variable name="node" select="."/>
  <!-- The next line is to deal with overloaded methods that have a const and a non-const version. -->
  <xsl:if test="not(@const=1) or empty($methods[@name=$node/@name and not(@const=1)])">
   <xsl:variable name="hide">
    <xsl:call-template name="should_hide"/>
   </xsl:variable>
   <xsl:if test="string-length($hide)=0">
    <xsl:call-template name="method">
     <xsl:with-param name="class_name" select="$class_name"/>
    </xsl:call-template>
   </xsl:if>
  </xsl:if>
 </xsl:for-each>
</xsl:template>

<!-- Print out the declaration for an enumeration -->
<xsl:template name="enumeration">
 <xsl:param name="class_name"/>
  <xsl:for-each select="EnumValue">
    integer*4, parameter :: OpenMM_<xsl:value-of select="concat($class_name, '_', @name, ' = ', @init)"/>
  </xsl:for-each>
 <xsl:value-of select="$newline"/>
</xsl:template>

<!-- Print out the declaration for a constructor -->
<xsl:template name="constructor">
 <xsl:param name="suffix"/>
        subroutine OpenMM_<xsl:value-of select="concat(@name, '_create', $suffix, '(result')"/>
  <xsl:for-each select="Argument">
   <xsl:value-of select="concat(', ', @name)"/>
  </xsl:for-each>
  <xsl:value-of select="')'"/>
            use OpenMM_Types; implicit none
            type (OpenMM_<xsl:value-of select="concat(@name, ') result', $newline)"/>
  <xsl:for-each select="Argument">
   <xsl:value-of select="'            '"/>
   <xsl:call-template name="declare_argument">
    <xsl:with-param name="type_id" select="@type"/>
    <xsl:with-param name="value" select="@name"/>
   </xsl:call-template>
   <xsl:value-of select="$newline"/>
  </xsl:for-each>
 <xsl:value-of select="'        end subroutine'"/>
</xsl:template>

<!-- Print out the declaration for a destructor -->
<xsl:template name="destructor">
        subroutine OpenMM_<xsl:value-of select="concat(@name, '_destroy(destroy)')"/>
            use OpenMM_Types; implicit none
            type (OpenMM_<xsl:value-of select="concat(@name, ') destroy', $newline)"/>
 <xsl:value-of select="'        end subroutine'"/>
</xsl:template>

<!-- Print out the declaration for a method -->
<xsl:template name="method">
 <xsl:param name="class_name"/>
 <xsl:variable name="has_return" select="@returns=$int_type_id or @returns=$double_type_id"/>
 <xsl:variable name="has_return_arg" select="not($has_return or @returns=$void_type_id)"/>
 <xsl:variable name="method" select="."/>
 <xsl:choose>
  <xsl:when test="$has_return">
   <xsl:text>
        function </xsl:text>
  </xsl:when>
  <xsl:otherwise>
   <xsl:text>
        subroutine </xsl:text>
  </xsl:otherwise>
 </xsl:choose>
 <xsl:value-of select="concat('OpenMM_', $class_name, '_', @name, '(')"/>
 <xsl:if test="not(@static='1')">
  <xsl:value-of select="'target'"/>
 </xsl:if>
 <xsl:for-each select="Argument">
  <xsl:if test="position() > 1 or not($method/@static='1')">
   <xsl:text>, &amp;
                          </xsl:text>
  </xsl:if>
  <xsl:value-of select="@name"/>
 </xsl:for-each>
 <xsl:if test="$has_return_arg">
  <xsl:if test="not(@static='1') or not(empty(Argument))">
   <xsl:text>, &amp;
                          </xsl:text>
  </xsl:if>
  <xsl:value-of select="'result'"/>
 </xsl:if>
 <xsl:value-of select="')'"/>
 <xsl:text>
            use OpenMM_Types; implicit none</xsl:text>
 <!-- Generate the list of argument types -->
 <xsl:if test="not(@static='1')">
            type (OpenMM_<xsl:value-of select="concat($class_name, ') target')"/>
 </xsl:if>
 <xsl:value-of select="$newline"/>
 <xsl:for-each select="Argument">
  <xsl:value-of select="'            '"/>
   <xsl:call-template name="declare_argument">
    <xsl:with-param name="type_id" select="@type"/>
    <xsl:with-param name="value" select="@name"/>
   </xsl:call-template>
   <xsl:value-of select="$newline"/>
 </xsl:for-each>
 <xsl:if test="$has_return">
  <xsl:value-of select="'            '"/>
   <xsl:call-template name="declare_argument">
    <xsl:with-param name="type_id" select="@returns"/>
    <xsl:with-param name="value" select="concat(' OpenMM_', $class_name, '_', @name)"/>
   </xsl:call-template>
  <xsl:value-of select="$newline"/>
 </xsl:if>
 <xsl:if test="$has_return_arg">
  <xsl:value-of select="'            '"/>
   <xsl:call-template name="declare_argument">
    <xsl:with-param name="type_id" select="@returns"/>
    <xsl:with-param name="value" select="'result'"/>
   </xsl:call-template>
  <xsl:value-of select="$newline"/>
 </xsl:if>
 <xsl:choose>
  <xsl:when test="$has_return">
   <xsl:value-of select="'        end function'"/>
  </xsl:when>
  <xsl:otherwise>
   <xsl:value-of select="'        end subroutine'"/>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

<!-- Print out the description of an argument -->
<xsl:template name="declare_argument">
 <xsl:param name="type_id"/>
 <xsl:param name="value"/>
 <xsl:variable name="node" select="/GCC_XML/*[@id=$type_id]"/>
 <xsl:choose>
  <xsl:when test="$type_id=$int_type_id">
   <xsl:value-of select="concat('integer*4 ', $value)"/>
  </xsl:when>
  <xsl:when test="$type_id=$bool_type_id">
   <xsl:value-of select="concat('integer*4 ', $value)"/>
  </xsl:when>
  <xsl:when test="$type_id=$double_type_id">
   <xsl:value-of select="concat('real*8 ', $value)"/>
  </xsl:when>
  <xsl:when test="$type_id=$string_type_id or $type_id=$ptr_const_char_type_id">
   <xsl:value-of select="concat('character(*) ', $value)"/>
  </xsl:when>
  <xsl:when test="$type_id=$vec3_type_id">
   <xsl:value-of select="concat('real*8 ', $value, '(3)')"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_string_type_id">
   <xsl:value-of select="concat('type (OpenMM_StringArray) ', $value)"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_vec3_type_id">
   <xsl:value-of select="concat('type (OpenMM_Vec3Array) ', $value)"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_bond_type_id">
   <xsl:value-of select="concat('type (OpenMM_BondArray) ', $value)"/>
  </xsl:when>
  <xsl:when test="$type_id=$map_parameter_type_id">
   <xsl:value-of select="concat('type (OpenMM_ParameterArray) ', $value)"/>
  </xsl:when>
  <xsl:when test="$type_id=$map_property_type_id">
   <xsl:value-of select="concat('type (OpenMM_PropertyArray) ', $value)"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_double_type_id">
   <xsl:value-of select="concat('type (OpenMM_DoubleArray) ', $value)"/>
  </xsl:when>
  <xsl:when test="$type_id=$vector_int_type_id">
   <xsl:value-of select="concat('type (OpenMM_IntArray) ', $value)"/>
  </xsl:when>
  <xsl:when test="local-name($node)='ReferenceType' or local-name($node)='PointerType'">
   <xsl:call-template name="declare_argument">
    <xsl:with-param name="type_id" select="$node/@type"/>
    <xsl:with-param name="value" select="$value"/>
   </xsl:call-template>
  </xsl:when>
  <xsl:when test="local-name($node)='CvQualifiedType'">
   <xsl:call-template name="declare_argument">
    <xsl:with-param name="type_id" select="$node/@type"/>
    <xsl:with-param name="value" select="$value"/>
   </xsl:call-template>
  </xsl:when>
  <xsl:when test="local-name($node)='Enumeration'">
   <xsl:value-of select="concat('integer*4 ', $value)"/>
  </xsl:when>
  <xsl:when test="$node/@context=$openmm_namespace_id">
   <xsl:value-of select="concat('type (OpenMM_', $node/@name, ') ', $value)"/>
  </xsl:when>
 </xsl:choose>
</xsl:template>

<!-- Determine whether a method should be hidden -->
<xsl:template name="should_hide">
 <xsl:variable name="class_id" select="@context"/>
 <xsl:variable name="method_name" select="concat('OpenMM_', /GCC_XML/Class[@id=$class_id]/@name, '_', @name)"/>
 <xsl:if test="not(empty(index-of($skip_methods, $method_name)))">1</xsl:if>
 <xsl:call-template name="hide_type">
  <xsl:with-param name="type_id" select="@returns"/>
 </xsl:call-template>
 <xsl:for-each select="Argument">
  <xsl:call-template name="hide_type">
   <xsl:with-param name="type_id" select="@type"/>
  </xsl:call-template>
 </xsl:for-each>
</xsl:template>

<!-- This is called by should_hide.  It generates output if the specified type should be hidden. -->
<xsl:template name="hide_type">
 <xsl:param name="type_id"/>
 <xsl:variable name="node" select="/GCC_XML/*[@id=$type_id]"/>
 <xsl:choose>
  <xsl:when test="local-name($node)='ReferenceType' or local-name($node)='PointerType' or local-name($node)='CvQualifiedType'">
   <xsl:call-template name="hide_type">
    <xsl:with-param name="type_id" select="$node/@type"/>
   </xsl:call-template>
  </xsl:when>
  <xsl:when test="local-name($node)='Enumeration'">
   <xsl:call-template name="hide_type">
    <xsl:with-param name="type_id" select="$node/@context"/>
   </xsl:call-template>
  </xsl:when>
  <xsl:when test="$node/@context=$openmm_namespace_id and not(empty(index-of($hide_classes, $node/@name)))">
   <xsl:value-of select="1"/>
  </xsl:when>
 </xsl:choose>
</xsl:template>

</xsl:stylesheet>
