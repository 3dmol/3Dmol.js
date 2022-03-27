(function (global, factory) {
  (factory((global['MMTF'] = global.MMTF || {})));
}(this, function (exports) { 'use strict';

  /**
   * @file utf8-utils
   * @private
   * @author Alexander Rose <alexander.rose@weirdbyte.de>
   * mostly copied from https://github.com/creationix/msgpack-js-browser
   * by Tim Caswell <tim@creationix.com>, MIT License, Copyright (c) 2013
   */


  // Encode string as utf8 into dataview at offset
  function utf8Write(view, offset, string) {
    var byteLength = view.byteLength;
    for(var i = 0, l = string.length; i < l; i++) {
      var codePoint = string.charCodeAt(i);

      // One byte of UTF-8
      if (codePoint < 0x80) {
        view.setUint8(offset++, codePoint >>> 0 & 0x7f | 0x00);
        continue;
      }

      // Two bytes of UTF-8
      if (codePoint < 0x800) {
        view.setUint8(offset++, codePoint >>> 6 & 0x1f | 0xc0);
        view.setUint8(offset++, codePoint >>> 0 & 0x3f | 0x80);
        continue;
      }

      // Three bytes of UTF-8.
      if (codePoint < 0x10000) {
        view.setUint8(offset++, codePoint >>> 12 & 0x0f | 0xe0);
        view.setUint8(offset++, codePoint >>> 6  & 0x3f | 0x80);
        view.setUint8(offset++, codePoint >>> 0  & 0x3f | 0x80);
        continue;
      }

      // Four bytes of UTF-8
      if (codePoint < 0x110000) {
        view.setUint8(offset++, codePoint >>> 18 & 0x07 | 0xf0);
        view.setUint8(offset++, codePoint >>> 12 & 0x3f | 0x80);
        view.setUint8(offset++, codePoint >>> 6  & 0x3f | 0x80);
        view.setUint8(offset++, codePoint >>> 0  & 0x3f | 0x80);
        continue;
      }
      throw new Error("bad codepoint " + codePoint);
    }
  }

  function utf8ByteCount(string) {
    var count = 0;
    for(var i = 0, l = string.length; i < l; i++) {
      var codePoint = string.charCodeAt(i);
      if (codePoint < 0x80) {
        count += 1;
        continue;
      }
      if (codePoint < 0x800) {
        count += 2;
        continue;
      }
      if (codePoint < 0x10000) {
        count += 3;
        continue;
      }
      if (codePoint < 0x110000) {
        count += 4;
        continue;
      }
      throw new Error("bad codepoint " + codePoint);
    }
    return count;
  }

  /**
   * encode data value (recursively) into binary encoded MessagePack v5 (http://msgpack.org/)
   * @param  {Object|Array|String|Number|Boolean|null} value  [description]
   * @param  {DataView} view   [description]
   * @param  {Integer} offset [description]
   * @return {Integer} number of bytes written into view
   */
  function encode$1(value, view, offset) {
    var type = typeof value;

    // Strings Bytes
    if (type === "string") {
      var length = utf8ByteCount(value);
      // fix str
      if (length < 0x20) {
        view.setUint8(offset, length | 0xa0);
        utf8Write(view, offset + 1, value);
        return 1 + length;
      }
      // str 8
      if (length < 0x100) {
        view.setUint8(offset, 0xd9);
        view.setUint8(offset + 1, length);
        utf8Write(view, offset + 2, value);
        return 2 + length;
      }
      // str 16
      if (length < 0x10000) {
        view.setUint8(offset, 0xda);
        view.setUint16(offset + 1, length);
        utf8Write(view, offset + 3, value);
        return 3 + length;
      }
      // str 32
      if (length < 0x100000000) {
        view.setUint8(offset, 0xdb);
        view.setUint32(offset + 1, length);
        utf8Write(view, offset + 5, value);
        return 5 + length;
      }
    }

    if (value instanceof Uint8Array) {
      var length = value.byteLength;
      var bytes = new Uint8Array(view.buffer);
      // bin 8
      if (length < 0x100) {
        view.setUint8(offset, 0xc4);
        view.setUint8(offset + 1, length);
        bytes.set(value, offset + 2);
        return 2 + length;
      }
      // bin 16
      if (length < 0x10000) {
        view.setUint8(offset, 0xc5);
        view.setUint16(offset + 1, length);
        bytes.set(value, offset + 3);
        return 3 + length;
      }
      // bin 32
      if (length < 0x100000000) {
        view.setUint8(offset, 0xc6);
        view.setUint32(offset + 1, length);
        bytes.set(value, offset + 5);
        return 5 + length;
      }
    }

    if (type === "number") {
      if (!isFinite(value)) {
          throw new Error("Number not finite: " + value);
      }

      // Floating point
      if (Math.floor(value) !== value) {
        view.setUint8(offset, 0xcb);
        view.setFloat64(offset + 1, value);
        return 9;
      }

      // Integers
      if (value >=0) {
        // positive fixnum
        if (value < 0x80) {
          view.setUint8(offset, value);
          return 1;
        }
        // uint 8
        if (value < 0x100) {
          view.setUint8(offset, 0xcc);
          view.setUint8(offset + 1, value);
          return 2;
        }
        // uint 16
        if (value < 0x10000) {
          view.setUint8(offset, 0xcd);
          view.setUint16(offset + 1, value);
          return 3;
        }
        // uint 32
        if (value < 0x100000000) {
          view.setUint8(offset, 0xce);
          view.setUint32(offset + 1, value);
          return 5;
        }
        throw new Error("Number too big 0x" + value.toString(16));
      }
      // negative fixnum
      if (value >= -0x20) {
        view.setInt8(offset, value);
        return 1;
      }
      // int 8
      if (value >= -0x80) {
        view.setUint8(offset, 0xd0);
        view.setInt8(offset + 1, value);
        return 2;
      }
      // int 16
      if (value >= -0x8000) {
        view.setUint8(offset, 0xd1);
        view.setInt16(offset + 1, value);
        return 3;
      }
      // int 32
      if (value >= -0x80000000) {
        view.setUint8(offset, 0xd2);
        view.setInt32(offset + 1, value);
        return 5;
      }
      throw new Error("Number too small -0x" + (-value).toString(16).substr(1));
    }

    // null
    if (value === null) {
      view.setUint8(offset, 0xc0);
      return 1;
    }

    // Boolean
    if (type === "boolean") {
      view.setUint8(offset, value ? 0xc3 : 0xc2);
      return 1;
    }

    // Container Types
    if (type === "object") {
      var length, size = 0;
      var isArray = Array.isArray(value);

      if (isArray) {
        length = value.length;
      }
      else {
        var keys = Object.keys(value);
        length = keys.length;
      }

      var size;
      if (length < 0x10) {
        view.setUint8(offset, length | (isArray ? 0x90 : 0x80));
        size = 1;
      }
      else if (length < 0x10000) {
        view.setUint8(offset, isArray ? 0xdc : 0xde);
        view.setUint16(offset + 1, length);
        size = 3;
      }
      else if (length < 0x100000000) {
        view.setUint8(offset, isArray ? 0xdd : 0xdf);
        view.setUint32(offset + 1, length);
        size = 5;
      }

      if (isArray) {
        for (var i = 0; i < length; i++) {
          size += encode$1(value[i], view, offset + size);
        }
      }
      else {
        for (var i = 0; i < length; i++) {
          var key = keys[i];
          size += encode$1(key, view, offset + size);
          size += encode$1(value[key], view, offset + size);
        }
      }

      return size;
    }
    throw new Error("Unknown type " + type);
  }

  function encodedSize(value) {
    var type = typeof value;

    // Raw Bytes
    if (type === "string") {
      var length = utf8ByteCount(value);
      if (length < 0x20) {
        return 1 + length;
      }
      if (length < 0x100) {
        return 2 + length;
      }
      if (length < 0x10000) {
        return 3 + length;
      }
      if (length < 0x100000000) {
        return 5 + length;
      }
    }

    if (value instanceof Uint8Array) {
      var length = value.byteLength;
      if (length < 0x100) {
        return 2 + length;
      }
      if (length < 0x10000) {
        return 3 + length;
      }
      if (length < 0x100000000) {
        return 5 + length;
      }
    }

    if (type === "number") {
      // Floating Point
      // double
      if (Math.floor(value) !== value) return 9;

      // Integers
      if (value >=0) {
        // positive fixnum
        if (value < 0x80) return 1;
        // uint 8
        if (value < 0x100) return 2;
        // uint 16
        if (value < 0x10000) return 3;
        // uint 32
        if (value < 0x100000000) return 5;
        throw new Error("Number too big 0x" + value.toString(16));
      }
      // negative fixnum
      if (value >= -0x20) return 1;
      // int 8
      if (value >= -0x80) return 2;
      // int 16
      if (value >= -0x8000) return 3;
      // int 32
      if (value >= -0x80000000) return 5;
      throw new Error("Number too small -0x" + value.toString(16).substr(1));
    }

    // Boolean, null
    if (type === "boolean" || value === null) return 1;

    // Container Types
    if (type === "object") {
      var length, size = 0;
      if (Array.isArray(value)) {
        length = value.length;
        for (var i = 0; i < length; i++) {
          size += encodedSize(value[i]);
        }
      }
      else {
        var keys = Object.keys(value);
        length = keys.length;
        for (var i = 0; i < length; i++) {
          var key = keys[i];
          size += encodedSize(key) + encodedSize(value[key]);
        }
      }
      if (length < 0x10) {
        return 1 + size;
      }
      if (length < 0x10000) {
        return 3 + size;
      }
      if (length < 0x100000000) {
        return 5 + size;
      }
      throw new Error("Array or object too long 0x" + length.toString(16));
    }
    throw new Error("Unknown type " + type);
  }

  function encodeMsgpack(value) {
    var buffer = new ArrayBuffer(encodedSize(value));
    var view = new DataView(buffer);
    encode$1(value, view, 0);
    return new Uint8Array(buffer);
  }

  /**
   * @file mmtf-constants
   * @private
   * @author Alexander Rose <alexander.rose@weirdbyte.de>
   */


  var PassThroughFields = [
      "mmtfVersion", "mmtfProducer",
      "unitCell", "spaceGroup", "structureId", "title",
      "depositionDate", "releaseDate",
      "experimentalMethods", "resolution", "rFree", "rWork",
      "bioAssemblyList", "ncsOperatorList", "entityList", "groupList",
      "numBonds", "numAtoms", "numGroups", "numChains", "numModels",
      "groupsPerChain", "chainsPerModel",
  ];

  var EncodedFields = [
  	// required
      "xCoordList", "yCoordList", "zCoordList",
      "groupIdList", "groupTypeList",
      "chainIdList",
      // optional
      "bFactorList", "atomIdList", "altLocList", "occupancyList",
      "secStructList", "insCodeList", "sequenceIndexList",
      "chainNameList",
      "bondAtomList", "bondOrderList"
  ];

  var AllFields = PassThroughFields.concat( EncodedFields );

  /**
   * @file mmtf-utils
   * @private
   * @author Alexander Rose <alexander.rose@weirdbyte.de>
   */

  /**
   * mmtf utils module.
   * @module MmtfUtils
   */


  function getView( ctor, typedArray, elemSize ){
      return typedArray ? new ctor(
          typedArray.buffer,
          typedArray.byteOffset,
          typedArray.byteLength / ( elemSize || 1 )
      ) : undefined;
  }

  function getDataView( typedArray ){
      return getView( DataView, typedArray );
  }

  /**
   * get an Uint8Array view on the input array memory
   * @static
   * @param  {TypedArray} dataArray - input array
   * @return {Uint8Array} new view on the input array memory
   */
  function getUint8View( typedArray ){
      return getView( Uint8Array, typedArray );
  }

  /**
   * get an Int8Array view on the input array memory
   * @static
   * @param  {TypedArray} dataArray - input array
   * @return {Int8Array} new view on the input array memory
   */
  function getInt8View( typedArray ){
      return getView( Int8Array, typedArray );
  }

  /**
   * get an Int32Array view on the input array memory
   * @static
   * @param  {TypedArray} dataArray - input array
   * @return {Int32Array} new view on the input array memory
   */
  function getInt32View( typedArray ){
      return getView( Int32Array, typedArray, 4 );
  }

  function getFloat32View( typedArray ){
      return getView( Float32Array, typedArray, 4 );
  }


  /**
   * get an Int16Array copy of the the input array data
   * @static
   * @param  {TypedArray} view - input data in big endian format
   * @param  {Int16Array} [dataArray] - pre-allocated output array
   * @return {Int16Array} copy of the input array data
   */
  function decodeInt16( bytes, output ){
      var n = bytes.length / 2;
      if( !output ) output = new Int16Array( n );
      for( var i = 0, i2 = 0; i < n; ++i, i2 += 2 ){
          output[ i ] = bytes[ i2 ] << 8 ^ bytes[ i2 + 1 ] << 0;
      }
      return output;
  }

  /**
   * make big endian buffer of an int16 array
   * @static
   * @param  {Array|TypedArray} array - array of int16 values
   * @return {ArrayBuffer} big endian buffer
   */
  function encodeInt16( array, output ){
      var n = array.length;
      if( !output ) output = new Uint8Array( 2 * n );
      var dv = getDataView( output );
      for( var i = 0; i < n; ++i ){
          dv.setInt16( 2 * i, array[ i ] );
      }
      return getUint8View( output );
  }

  /**
   * get an Int32Array copy of the the input array data
   * @static
   * @param  {TypedArray} view - input data in big endian format
   * @param  {Int32Array} [dataArray] - pre-allocated output array
   * @return {Int32Array} copy of the input array data
   */
  function decodeInt32( bytes, output ){
      var n = bytes.length / 4;
      if( !output ) output = new Int32Array( n );
      for( var i = 0, i4 = 0; i < n; ++i, i4 += 4 ){
          output[ i ] = (
              bytes[ i4     ] << 24 ^ bytes[ i4 + 1 ] << 16 ^
              bytes[ i4 + 2 ] <<  8 ^ bytes[ i4 + 3 ] <<  0
          );
      }
      return output;
  }

  /**
   * make big endian buffer of an int32 array
   * @static
   * @param  {Array|TypedArray} array - array of int32 values
   * @return {ArrayBuffer} big endian buffer
   */
  function encodeInt32( array, output ){
      var n = array.length;
      if( !output ) output = new Uint8Array( 4 * n );
      var dv = getDataView( output );
      for( var i = 0; i < n; ++i ){
          dv.setInt32( 4 * i, array[ i ] );
      }
      return getUint8View( output );
  }

  function decodeFloat32( bytes, output ){
      var n = bytes.length;
      if( !output ) output = new Float32Array( n / 4 );
      var dvOut = getDataView( output );
      var dvIn = getDataView( bytes );
      for( var i = 0, i4 = 0, il = n / 4; i < il; ++i, i4 += 4 ){
          dvOut.setFloat32( i4, dvIn.getFloat32( i4 ), true );
      }
      return output;
  }

  /**
   * decode integers into floats using given divisor
   * example:
   *     intArray: [ 12, 34, 543, 687, 2, 0, 4689 ]
   *     divisor: 100
   *     return: [ 0.12, 0.34, 5.43, 6.87, 0.02, 0.00, 46.89 ]
   * @static
   * @param  {TypedArray|Array} intArray - input array containing integers
   * @param  {Number} divisor - number to devide the integers to obtain floats
   * @param  {Float32Array} [dataArray] - pre-allocated output array
   * @return {Float32Array} decoded array
   */
  function decodeInteger( intArray, divisor, output ){
      var n = intArray.length;
      var invDiv = 1/divisor;
      if( !output ) output = new Float32Array( n );
      for( var i = 0; i < n; ++i ){
          // multiply by inverse of the divisor which is faster then division
          output[ i ] = intArray[ i ] * invDiv;
      }
      return output;
  }

  function encodeInteger( floatArray, factor, output ){
      var n = floatArray.length;
      if( !output ) output = new Int32Array( n );
      for( var i = 0; i < n; ++i ){
          output[ i ] = Math.round( floatArray[ i ] * factor );
      }
      return output;
  }



  /**
   * perform run-length decoding of input array
   * example:
   *     array: [ 0, 2, 3, 5 ]  // pairs of values and length of a run
   *     return: [ 0, 0, 3, 3, 3, 3, 3 ]
   * @static
   * @param  {TypedArray|Array} array - run-length encoded input array
   * @param  {TypedArray|Array} [dataArray] - pre-allocated output array
   * @return {TypedArray|Array} decoded array
   */
  function decodeRun( array, output ){
      var i, il;
      if( !output ){
          // calculate the length the decoded array will have
          var fullLength = 0;
          for( i = 0, il = array.length; i < il; i+=2 ){
              fullLength += array[ i + 1 ];
          }
          // create a new array of the same type of the input array
          output = new array.constructor( fullLength );
      }
      var dataOffset = 0;
      for( i = 0, il = array.length; i < il; i+=2 ){
          var value = array[ i ];  // value to be repeated
          var length = array[ i + 1 ];  // number of repeats
          for( var j = 0; j < length; ++j ){
              output[ dataOffset ] = value;
              ++dataOffset;
          }
      }
      return output;
  }

  function encodeRun( array ){
      if( array.length === 0 ) return new Int32Array();
      var i, il;
      // calculate output size
      var fullLength = 2;
      for( i = 1, il = array.length; i < il; ++i ){
          if( array[ i - 1 ] !== array[ i ] ){
              fullLength += 2;
          }
      }
      var output = new Int32Array( fullLength );
      var offset = 0;
      var runLength = 1;
      for( i = 1, il = array.length; i < il; ++i ){
          if( array[ i - 1 ] !== array[ i ] ){
              output[ offset ] = array[ i - 1 ];
              output[ offset + 1 ] = runLength;
              runLength = 1;
              offset += 2;
          }else{
              ++runLength;
          }
      }
      output[ offset ] = array[ array.length - 1 ];
      output[ offset + 1 ] = runLength;
      return output;
  }



  /**
   * perform delta decoding of the input array
   * by iterativly adding the ith element's value to the i+1th
   * example:
   *     dataArray: [ 0, 2, 1, 2, 1, 1, -4, -2, 9 ]
   *     return: [ 0, 2, 3, 5, 6, 7, 3, 1, 10 ]
   * @static
   * @param  {TypedArray|Array} dataArray - delta encoded input array
   * @return {TypedArray|Array} decoded array
   */
  function decodeDelta( array, output ){
      var n = array.length;
      if( !output ) output = new array.constructor( n );
      if( n ) output[ 0 ] = array[ 0 ];
      for( var i = 1; i < n; ++i ){
          output[ i ] = array[ i ] + output[ i - 1 ];
      }
      return output;
  }

  function encodeDelta( array, output ){
      var n = array.length;
      if( !output ) output = new array.constructor( n );
      output[ 0 ] = array[ 0 ];
      for( var i = 1; i < n; ++i ){
          output[ i ] = array[ i ] - array[ i - 1 ];
      }
      return output;
  }



  /**
   * [decodePacking description]
   * @param  {Int16Array|Int8Array} int16or8 [description]
   * @param  {Int32Array} output   [description]
   * @return {Int32Array}          [description]
   */
  function decodePacking( int16or8, output ){
      var upperLimit = int16or8 instanceof Int8Array ? 0x7F : 0x7FFF;
      var lowerLimit = -upperLimit - 1;
      var n = int16or8.length;
      var i, j;
      if( !output ){
          var fullLength = 0;
          for( i = 0; i < n; ++i ){
              if( int16or8[ i ] < upperLimit && int16or8[ i ] > lowerLimit ){
                  ++fullLength;
              }
          }
          output = new Int32Array( fullLength );
      }
      i = 0;
      j = 0;
      while( i < n ){
          var value = 0;
          while( int16or8[ i ] === upperLimit || int16or8[ i ] === lowerLimit ){
              value += int16or8[ i ];
              ++i;
          }
          value += int16or8[ i ];
          ++i;
          output[ j ] = value;
          ++j;
      }
      return output;
  }

  /**
   * integer packing using recursive indexing
   * @param  {Array|TyepedArray} intArray [description]
   * @param  {Boolean} useInt8  [description]
   * @return {Int16Array|Int8Array}          [description]
   */
  function encodePacking( intArray, useInt8 ){
      var upperLimit = useInt8 ? 0x7F : 0x7FFF;
      var lowerLimit = -upperLimit - 1;
      var i;
      var n = intArray.length;
      var size = 0;
      for( i = 0; i < n; ++i ){
          var value = intArray[ i ];
          if( value === 0 ){
              ++size;
          }else if( value === upperLimit || value === lowerLimit ){
              size += 2;
          }else if( value > 0) {
              size += Math.ceil( value / upperLimit );
          }else {
              size += Math.ceil( value / lowerLimit );
          }
      }
      var output = useInt8 ? new Int8Array( size ) : new Int16Array( size );
      var j = 0;
      for( i = 0; i < n; ++i ){
          var value = intArray[ i ];
          if( value >= 0) {
              while( value >= upperLimit ){
                  output[ j ] = upperLimit;
                  ++j;
                  value -= upperLimit;
              }
          }else{
              while( value <= lowerLimit ){
                  output[ j ] = lowerLimit;
                  ++j;
                  value -= lowerLimit;
              }
          }
          output[ j ] = value;
          ++j;
      }
      return output;
  }



  function decodeDeltaRun( array, output ){
      return decodeDelta( decodeRun( array ), output );
  }

  function encodeDeltaRun( array ){
      return encodeRun( encodeDelta( array ) );
  }



  /**
   * perform run-length decoding followed (@see decodeRunLength)
   * by decoding integers into floats using given divisor (@see decodeIntegerToFloat)
   * example:
   *     array: [ 320, 3, 100, 2 ]
   *     divisor: 100
   *     return: [ 3.20, 3.20, 3.20, 1.00, 1.00 ]
   * @static
   * @param  {Uint8Array} array - run-length encoded int32 array as bytes in big endian format
   * @param  {Integer} divisor - number to devide the integers to obtain floats
   * @param  {Float32Array} dataArray - pre-allocated output array
   * @return {Float32Array} decoded array
   */
  function decodeIntegerRun( intArray, divisor, output ){
      return decodeInteger(
          decodeRun( intArray, getInt32View( output ) ), divisor, output
      );
  }

  function encodeIntegerRun( floatArray, factor ){
      return encodeRun( encodeInteger( floatArray, factor ) );
  }



  function decodeIntegerDelta( intArray, divisor, output ){
      return decodeInteger(
          decodeDelta( intArray, getInt32View( output ) ), divisor, output
      );
  }

  function encodeIntegerDelta( floatArray, factor, output ){
      return encodeDelta( encodeInteger( floatArray, factor ), output );
  }



  function decodeIntegerPacking( int16or8, divisor, output ){
      return decodeInteger(
          decodePacking( int16or8, getInt32View( output ) ), divisor, output
      );
  }

  function decodeIntegerDeltaPacking( int16or8, divisor, output ){
      var unpacked = decodePacking( int16or8, getInt32View( output ) );
      return decodeIntegerDelta( unpacked, divisor, getFloat32View( unpacked ) );
  }

  function encodeIntegerDeltaPacking( floatArray, factor, useInt8 ){
      return encodePacking( encodeIntegerDelta( floatArray, factor ), useInt8 );
  }



  function decodeBytes( bytes ){
      var dv = getDataView( bytes );
      var type = dv.getInt32( 0 );
      var size = dv.getInt32( 4 );
      var param = bytes.subarray( 8, 12 );
      var bytes = bytes.subarray( 12 );
      return [ type, bytes, size, param ];
  }

  function encodeBytes( type, size, param, bytes ){
      var buffer = new ArrayBuffer( 12 + bytes.byteLength );
      var out = new Uint8Array( buffer );
      var dv = new DataView( buffer );
      dv.setInt32( 0, type );
      dv.setInt32( 4, size );
      if( param ) out.set( param, 8 );
      out.set( bytes, 12 );
      return out;
  }

  function passInt8( int8 ){
      var size = int8.length;
      var bytes = getUint8View( int8 );
      return encodeBytes( 2, size, undefined, bytes );
  }

  function passInt32( int32 ){
      var size = int32.length;
      var bytes = encodeInt32( int32 );
      return encodeBytes( 4, size, undefined, bytes );
  }

  function passString( stringBytes, length ){
      var size = stringBytes.length / length;
      var param = encodeInt32([ length ]);
      var bytes = getUint8View( stringBytes );
      return encodeBytes( 5, size, param, bytes );
  }

  function runChar( charBytes ){
      var size = charBytes.length;
      var bytes = encodeInt32( encodeRun( charBytes ) );
      return encodeBytes( 6, size, undefined, bytes );
  }

  function deltaRun( int32 ){
      var size = int32.length;
      var bytes = encodeInt32( encodeDeltaRun( int32 ) );
      return encodeBytes( 8, size, undefined, bytes );
  }

  function integerRun( float32, factor ){
      var size = float32.length;
      var param = encodeInt32([ factor ]);
      var bytes = encodeInt32( encodeIntegerRun( float32, factor ) );
      return encodeBytes( 9, size, param, bytes );
  }

  function integerDeltaPacking16( float32, factor ){
      var size = float32.length;
      var param = encodeInt32([ factor ]);
      var bytes = encodeInt16( encodeIntegerDeltaPacking( float32, factor ) );
      return encodeBytes( 10, size, param, bytes );
  }

  function encodeMmtf( inputDict ){

      var outputDict = {};

      // copy some fields over from the input dict
      PassThroughFields.forEach( function( name ){
          if( inputDict[ name ] !== undefined ){
              outputDict[ name ] = inputDict[ name ];
          }
      } );

      //////////////
      // bond data

      // encode inter group bond atom indices, i.e. get bytes in big endian order
      if( inputDict.bondAtomList ){
          outputDict.bondAtomList = passInt32( inputDict.bondAtomList );
      }

      // encode inter group bond orders, i.e. get bytes
      if( inputDict.bondOrderList ){
          outputDict.bondOrderList = passInt8( inputDict.bondOrderList );
      }

      //////////////
      // atom data

      // split-list delta & integer encode x, y, z atom coords
      outputDict.xCoordList = integerDeltaPacking16( inputDict.xCoordList, 1000 );
      outputDict.yCoordList = integerDeltaPacking16( inputDict.yCoordList, 1000 );
      outputDict.zCoordList = integerDeltaPacking16( inputDict.zCoordList, 1000 );

      // split-list delta & integer encode b-factors
      if( inputDict.bFactorList ){
          outputDict.bFactorList = integerDeltaPacking16( inputDict.bFactorList, 100 );
      }

      // delta & run-length encode atom ids
      if( inputDict.atomIdList ){
          outputDict.atomIdList = deltaRun( inputDict.atomIdList );
      }

      // run-length encode alternate labels
      if( inputDict.altLocList ){
          outputDict.altLocList = runChar( inputDict.altLocList );
      }

      // run-length & integer encode occupancies
      if( inputDict.occupancyList ){
          outputDict.occupancyList = integerRun( inputDict.occupancyList, 100 );
      }

      ///////////////
      // group data

      // run-length & delta encode group numbers
      outputDict.groupIdList = deltaRun( inputDict.groupIdList );

      // encode group types, i.e. get int32 array
      outputDict.groupTypeList = passInt32( inputDict.groupTypeList );

      // encode secondary structure, i.e. get bytes
      if( inputDict.secStructList ){
          outputDict.secStructList = passInt8( inputDict.secStructList );
      }

      // run-length encode insertion codes
      if( inputDict.insCodeList ){
          outputDict.insCodeList = runChar( inputDict.insCodeList );
      }

      // run-length & delta encode sequence indices
      if( inputDict.sequenceIndexList ){
          outputDict.sequenceIndexList = deltaRun( inputDict.sequenceIndexList );
      }

      ///////////////
      // chain data

      // encode chain ids, i.e. get bytes
      outputDict.chainIdList = passString( inputDict.chainIdList, 4 );

      // encode chain names, i.e. get bytes
      if( inputDict.chainNameList ){
          outputDict.chainNameList = passString( inputDict.chainNameList, 4 );
      }

      return outputDict;

  }

  /**
   * @file msgpack-decode
   * @private
   * @author Alexander Rose <alexander.rose@weirdbyte.de>
   */

  /**
   * msgpack decode module.
   * @module MsgpackDecode
   */

  /**
   * decode binary encoded MessagePack v5 (http://msgpack.org/) data
   * @static
   * @param  {Uint8Array} buffer - binary encoded MessagePack data
   * @return {Object|Array|String|Number|Boolean|null} decoded Messagepack data
   */
  function decodeMsgpack(buffer) {
    // Loosely based on
    // The MIT License (MIT)
    // Copyright (c) 2013 Tim Caswell <tim@creationix.com>
    // https://github.com/creationix/msgpack-js
    var offset = 0;
    var dataView = new DataView(buffer.buffer);

    /**
     * decode all key-value pairs of a map into an object
     * @param  {Integer} length - number of key-value pairs
     * @return {Object} decoded map
     */
    function map(length) {
      var value = {};
      for (var i = 0; i < length; i++) {
        var key = parse();
        value[key] = parse();
      }
      return value;
    }

    /**
     * decode binary array
     * @param  {Integer} length - number of elements in the array
     * @return {Uint8Array} decoded array
     */
    function bin(length) {
      var value = buffer.subarray(offset, offset + length);
      offset += length;
      return value;
    }

    /**
     * decode string
     * @param  {Integer} length - number string characters
     * @return {String} decoded string
     */
    function str(length) {
      var array = buffer.subarray(offset, offset + length);
      offset += length;
      // limit number of arguments to String.fromCharCode to something
      // browsers can handle, see http://stackoverflow.com/a/22747272
      var chunkSize = 0xffff;
      if(length > chunkSize){
        var c = [];
        for(var i = 0; i < array.length; i += chunkSize) {
          c.push(String.fromCharCode.apply(
            null, array.subarray(i, i + chunkSize)
          ));
        }
        return c.join("");
      }else{
        return String.fromCharCode.apply(null, array);
      }
    }

    /**
     * decode array
     * @param  {Integer} length - number of array elements
     * @return {Array} decoded array
     */
    function array(length) {
      var value = new Array(length);
      for (var i = 0; i < length; i++) {
        value[i] = parse();
      }
      return value;
    }

    /**
     * recursively parse the MessagePack data
     * @return {Object|Array|String|Number|Boolean|null} decoded MessagePack data
     */
    function parse() {
      var type = buffer[offset];
      var value, length, extType;
      // Positive FixInt
      if ((type & 0x80) === 0x00) {
        offset++;
        return type;
      }
      // FixMap
      if ((type & 0xf0) === 0x80) {
        length = type & 0x0f;
        offset++;
        return map(length);
      }
      // FixArray
      if ((type & 0xf0) === 0x90) {
        length = type & 0x0f;
        offset++;
        return array(length);
      }
      // FixStr
      if ((type & 0xe0) === 0xa0) {
        length = type & 0x1f;
        offset++;
        return str(length);
      }
      // Negative FixInt
      if ((type & 0xe0) === 0xe0) {
        value = dataView.getInt8(offset);
        offset++;
        return value;
      }
      switch (type) {
      // nil
      case 0xc0:
        offset++;
        return null;
      // 0xc1: (never used, could be employed for padding)
      // false
      case 0xc2:
        offset++;
        return false;
      // true
      case 0xc3:
        offset++;
        return true;
      // bin 8
      case 0xc4:
        length = dataView.getUint8(offset + 1);
        offset += 2;
        return bin(length);
      // bin 16
      case 0xc5:
        length = dataView.getUint16(offset + 1);
        offset += 3;
        return bin(length);
      // bin 32
      case 0xc6:
        length = dataView.getUint32(offset + 1);
        offset += 5;
        return bin(length);
      // // ext 8
      // case 0xc7:
      //   length = dataView.getUint8(offset + 1);
      //   extType = dataView.getUint8(offset + 2);
      //   offset += 3;
      //   return [extType, bin(length)];
      // // ext 16
      // case 0xc8:
      //   length = dataView.getUint16(offset + 1);
      //   extType = dataView.getUint8(offset + 3);
      //   offset += 4;
      //   return [extType, bin(length)];
      // // ext 32
      // case 0xc9:
      //   length = dataView.getUint32(offset + 1);
      //   extType = dataView.getUint8(offset + 5);
      //   offset += 6;
      //   return [extType, bin(length)];
      // float 32
      case 0xca:
        value = dataView.getFloat32(offset + 1);
        offset += 5;
        return value;
      // float 64
      case 0xcb:
        value = dataView.getFloat64(offset + 1);
        offset += 9;
        return value;
      // uint8
      case 0xcc:
        value = buffer[offset + 1];
        offset += 2;
        return value;
      // uint 16
      case 0xcd:
        value = dataView.getUint16(offset + 1);
        offset += 3;
        return value;
      // uint 32
      case 0xce:
        value = dataView.getUint32(offset + 1);
        offset += 5;
        return value;
      // // uint64
      // case 0xcf:
      //   // FIXME not available/representable in JS
      //   // largest possible int in JS is 2^53
      //   // value = dataView.getUint64(offset + 1);
      //   offset += 9;
      //   return 0;
      // int 8
      case 0xd0:
        value = dataView.getInt8(offset + 1);
        offset += 2;
        return value;
      // int 16
      case 0xd1:
        value = dataView.getInt16(offset + 1);
        offset += 3;
        return value;
      // int 32
      case 0xd2:
        value = dataView.getInt32(offset + 1);
        offset += 5;
        return value;
      // // int 64
      // case 0xd3:
      //   // FIXME not available/representable in JS
      //   // largest possible int in JS is 2^53
      //   // value = dataView.getInt64(offset + 1);
      //   offset += 9;
      //   return 0;

      // // fixext 1
      // case 0xd4:
      //   extType = dataView.getUint8(offset + 1);
      //   offset += 2;
      //   return [extType, bin(1)];
      // // fixext 2
      // case 0xd5:
      //   extType = dataView.getUint8(offset + 1);
      //   offset += 2;
      //   return [extType, bin(2)];
      // // fixext 4
      // case 0xd6:
      //   extType = dataView.getUint8(offset + 1);
      //   offset += 2;
      //   return [extType, bin(4)];
      // // fixext 8
      // case 0xd7:
      //   extType = dataView.getUint8(offset + 1);
      //   offset += 2;
      //   return [extType, bin(8)];
      // // fixext 16
      // case 0xd8:
      //   extType = dataView.getUint8(offset + 1);
      //   offset += 2;
      //   return [extType, bin(16)];
      // str 8
      case 0xd9:
        length = dataView.getUint8(offset + 1);
        offset += 2;
        return str(length);
      // str 16
      case 0xda:
        length = dataView.getUint16(offset + 1);
        offset += 3;
        return str(length);
      // str 32
      case 0xdb:
        length = dataView.getUint32(offset + 1);
        offset += 5;
        return str(length);
      // array 16
      case 0xdc:
        length = dataView.getUint16(offset + 1);
        offset += 3;
        return array(length);
      // array 32
      case 0xdd:
        length = dataView.getUint32(offset + 1);
        offset += 5;
        return array(length);
      // map 16:
      case 0xde:
        length = dataView.getUint16(offset + 1);
        offset += 3;
        return map(length);
      // map 32
      case 0xdf:
        length = dataView.getUint32(offset + 1);
        offset += 5;
        return map(length);
      }

      throw new Error("Unknown type 0x" + type.toString(16));
    }

    // start the recursive parsing
    return parse();
  }

  /**
   * Fields shared in encoded and decoded mmtf data objects.
   * @typedef {Object} module:MmtfDecode.SharedMmtfData
   * @property {String} mmtfVersion - MMTF specification version
   * @property {String} mmtfProducer - Program that created the file
   * @property {Float[]} [unitCell] - Crystallographic unit cell
   * @property {Float} unitCell.0 - x length
   * @property {Float} unitCell.1 - y length
   * @property {Float} unitCell.2 - z length
   * @property {Float} unitCell.3 - alpha angle
   * @property {Float} unitCell.4 - beta angle
   * @property {Float} unitCell.5 - gamma angle
   * @property {String} [spaceGroup] - Hermann-Mauguin symbol
   * @property {String} [structureId] - Some reference, e.g. a PDB ID
   * @property {String} [title] - Short description
   * @property {String} [depositionDate] - Deposition date in YYYY-MM-DD format
   * @property {String} [releaseDate] - Release date in YYYY-MM-DD format
   * @property {String[]} [experimentalMethods] - Structure determination methods
   * @property {Float} [resolution] - Resolution in Å
   * @property {Float} [rFree] - R-free value
   * @property {Float} [rWork] - R-work value
   * @property {Integer} numBonds - Number of bonds
   * @property {Integer} numAtoms - Number of atoms
   * @property {Integer} numGroups - Number of groups (residues)
   * @property {Integer} numChains - Number of chains
   * @property {Integer} numModels - Number of models
   * @property {Integer[]} chainsPerModel - List of number of chains in each model
   * @property {Integer[]} groupsPerChain - List of number of groups in each chain
   * @property {Entity[]} [entityList] - List of entity objects
   * @property {Integer[]} entityList.chainIndexList - Pointers into chain data fields
   * @property {String} entityList.description - Description of the entity
   * @property {String} entityList.type - Name of the entity type
   * @property {String} entityList.sequence - One letter code sequence
   * @property {Assembly[]} [bioAssemblyList] - List of assembly objects
   * @property {Transform[]} bioAssemblyList.transformList - List of transform objects
   * @property {Integer[]} bioAssemblyList.transformList.chainIndexList - Pointers into chain data fields
   * @property {Float[]} bioAssemblyList.transformList.matrix - 4x4 transformation matrix
   * @property {Array[]} [ncsOperatorList] - List of ncs operator matrices
   * @property {Float[]} ncsOperatorList. - 4x4 transformation matrix
   * @property {GroupType[]} groupList - List of groupType objects
   * @property {Integer[]} groupList.formalChargeList - List of atom formal charges
   * @property {String[]} groupList.elementList - List of elements
   * @property {String[]} groupList.atomNameList - List of atom names
   * @property {Integer[]} groupList.bondAtomList - List of bonded atom indices
   * @property {Integer[]} groupList.bondOrderList - List of bond orders
   * @property {String} groupList.groupName - The name of the group
   * @property {String} groupList.singleLetterCode - The single letter code
   * @property {String} groupList.chemCompType -  The chemical component type
   */

  /**
   * Encoded mmtf data object. Also includes the fields from {@link module:MmtfDecode.SharedMmtfData}. See MMTF specification on how they are encoded.
   * @typedef {Object} module:MmtfDecode.EncodedMmtfData
   * @mixes module:MmtfDecode.SharedMmtfData
   * @property {Uint8Array} [bondAtomList] - Encoded bonded atom indices
   * @property {Uint8Array} [bondOrderList] - Encoded bond orders
   * @property {Uint8Array} xCoordBig - Encoded x coordinates in Å, part 1
   * @property {Uint8Array} xCoordSmall - Encoded x coordinates in Å, part 2
   * @property {Uint8Array} yCoordBig - Encoded y coordinates in Å, part 1
   * @property {Uint8Array} yCoordSmall - Encoded y coordinates in Å, part 2
   * @property {Uint8Array} yCoordBig - Encoded y coordinates in Å, part 1
   * @property {Uint8Array} yCoordSmall - Encoded y coordinates in Å, part 2
   * @property {Uint8Array} [bFactorBig] - Encoded B-factors in Å^2, part 1
   * @property {Uint8Array} [bFactorSmall] - Encoded B-factors in Å^2, part 2
   * @property {Uint8Array} [atomIdList] - Encoded  atom ids
   * @property {Uint8Array} [altLocList] - Encoded alternate location labels
   * @property {Uint8Array} [occupancyList] - Encoded occupancies
   * @property {Uint8Array} groupIdList - Encoded group ids
   * @property {Uint8Array} groupTypeList - Encoded group types
   * @property {Uint8Array} [secStructList] - Encoded secondary structure codes
   * @property {Uint8Array} [insCodeList] - Encoded insertion codes
   * @property {Uint8Array} [seuenceIdList] - Encoded sequence ids
   * @property {Uint8Array} chainIdList - Encoded chain ids
   * @property {Uint8Array} [chainNameList] - Encoded chain names
   */

  /**
   * Decoded mmtf data object. Also includes fields the from {@link module:MmtfDecode.SharedMmtfData}.
   * @typedef {Object} module:MmtfDecode.MmtfData
   * @mixes module:MmtfDecode.SharedMmtfData
   * @property {Int32Array} [bondAtomList] - List of bonded atom indices
   * @property {Uint8Array} [bondOrderList] - List of bond orders
   * @property {Float32Array} xCoordList - List of x coordinates in Å
   * @property {Float32Array} yCoordList - List of y coordinates in Å
   * @property {Float32Array} zCoordList - List of z coordinates in Å
   * @property {Float32Array} [bFactorList] - List of B-factors in Å^2
   * @property {Int32Array} [atomIdList] - List of atom ids
   * @property {Uint8Array} [altLocList] - List of alternate location labels
   * @property {Float32Array} [occupancyList] - List of occupancies
   * @property {Int32Array} groupIdList - List of group ids
   * @property {Int32Array} groupTypeList - List of group types
   * @property {Int8Array} [secStructList] - List of secondary structure codes, encoding
   *    0: pi helix, 1: bend, 2: alpha helix, 3: extended,
   *    4: 3-10 helix, 5: bridge, 6: turn, 7: coil, -1: undefined
   * @property {Uint8Array} [insCodeList] - List of insertion codes
   * @property {Int32Array} [seuenceIdList] - List of sequence ids
   * @property {Uint8Array} chainIdList - List of chain ids
   * @property {Uint8Array} [chainNameList] - List of chain names
   */


  /**
   * [performDecoding description]
   * @param  {Integer} bytes [description]
   * @param  {Integer} size  [description]
   * @param  {Uint8Array} param [description]
   * @return {TypedArray}       [description]
   */
  function performDecoding( type, bytes, size, param ){

      switch( type ){
          case 1:
              return decodeFloat32( bytes );
          case 2:
              return getInt8View( bytes );
          case 3:
              return decodeInt16( bytes );
          case 4:
              return decodeInt32( bytes );
          case 5:
              // var length = decodeInt32( param )[ 0 ];
              return getUint8View( bytes );  // interpret as string array
          case 6:
              // interpret as char array
              return decodeRun( decodeInt32( bytes ), new Uint8Array( size ) );
          case 7:
              return decodeRun( decodeInt32( bytes ) )
          case 8:
              return decodeDeltaRun( decodeInt32( bytes ) );
          case 9:
              return decodeIntegerRun( decodeInt32( bytes ), decodeInt32( param )[ 0 ] );
          case 10:
              return decodeIntegerDeltaPacking( decodeInt16( bytes ), decodeInt32( param )[ 0 ] );
          case 11:
              return decodeInteger( decodeInt16( bytes ), decodeInt32( param )[ 0 ] );
          case 12:
              return decodeIntegerPacking( decodeInt16( bytes ), decodeInt32( param )[ 0 ] );
          case 13:
              return decodeIntegerPacking( getInt8View( bytes ), decodeInt32( param )[ 0 ] );
          case 14:
              return decodePacking( decodeInt16( bytes ) );
          case 15:
              return decodePacking( getInt8View( bytes ) );
      }

  };


  /**
   * Decode MMTF fields
   * @static
   * @param  {Object} inputDict - encoded MMTF data
   * @param  {Object} [params] - decoding parameters
   * @param  {String[]} params.ignoreFields - names of optional fields not to decode
   * @return {module:MmtfDecode.MmtfData} mmtfData
   */
  function decodeMmtf( inputDict, params ){

      params = params || {};
      var ignoreFields = params.ignoreFields;
      var outputDict = {};

      AllFields.forEach( function( name ){
          var ignore = ignoreFields ? ignoreFields.indexOf( name ) !== -1 : false;
          var data = inputDict[ name ];
          if( !ignore && data !== undefined ){
              if( data instanceof Uint8Array ){
                  outputDict[ name ] = performDecoding.apply( null, decodeBytes( data ) );
              }else{
                  outputDict[ name ] = data;
              }
          }
      } );

      return outputDict;

  }

  /**
   * @file mmtf-traverse
   * @private
   * @author Alexander Rose <alexander.rose@weirdbyte.de>
   */

  /**
   * mmtf traverse module.
   * @module MmtfTraverse
   */

  /**
   * Converts an array of ASCII codes trimming '\0' bytes
   * @private
   * @param  {Array} charCodeArray - array of ASCII char codes
   * @return {String} '\0' trimmed string
   */
  function fromCharCode( charCodeArray ){
      return String.fromCharCode.apply( null, charCodeArray ).replace(/\0/g, '');
  }


  /**
   * @callback module:MmtfTraverse.onModel
   * @param {Object} modelData
   * @param {Integer} modelData.chainCount - number of chains in the model
   * @param {Integer} modelData.modelIndex - index of the model
   */

  /**
   * @callback module:MmtfTraverse.onChain
   * @param {Object} chainData
   * @param {Integer} chainData.groupCount - number of groups in the chain
   * @param {Integer} chainData.chainIndex - index of the chain
   * @param {Integer} chainData.modelIndex - index of the parent model
   * @param {String} chainData.chainId - chain id
   * @param {?String} chainData.chainName - additional chain name
   */

  /**
   * @callback module:MmtfTraverse.onGroup
   * @param {Object} groupData
   * @param {Integer} groupData.atomCount - number of atoms in the group
   * @param {Integer} groupData.groupIndex - index of the group
   * @param {Integer} groupData.chainIndex - index of the parent chain
   * @param {Integer} groupData.modelIndex - index of the parent model
   * @param {Integer} groupData.groupId - group id (residue number)
   * @param {Integer} groupData.groupType - index to an entry in {@link module:MmtfDecode.MmtfData}#groupList
   * @param {String} groupData.groupName - name of the group, 0 to 5 characters
   * @param {Char} groupData.singleLetterCode - IUPAC single letter code, otherwise 'X', 1 character
   * @param {String} groupData.chemCompType - chemical component type from the mmCIF dictionary
   * @param {?Integer} groupData.secStruct - sencoded secondary structure |
   *    0: pi helix, 1: bend, 2: alpha helix, 3: extended,
   *    4: 3-10 helix, 5: bridge, 6: turn, 7: coil, -1: undefined
   * @param {?Char} groupData.insCode - insertion code
   * @param {?Integer} groupData.sequenceIndex - index to the `sequence` property of
   *    the corresponding entity, -1 when the entity has no sequence
   */

  /**
   * @callback module:MmtfTraverse.onAtom
   * @param {Object} atomData
   * @param {Integer} atomData.atomIndex - index of the atom
   * @param {Integer} atomData.groupIndex - index of the parent group
   * @param {Integer} atomData.chainIndex - index of the parent chain
   * @param {Integer} atomData.modelIndex - index of the parent model
   * @param {?Integer} atomData.atomId - atom id
   * @param {String} atomData.element - IUPAC element name, 0 to 3 characters
   * @param {String} atomData.atomName - name of the atom, 0 to 5 characters
   * @param {Integer} atomData.formalCharge - formal charge of the atom
   * @param {Float} atomData.xCoord - x coordinate in Å
   * @param {Float} atomData.yCoord - y coordinate in Å
   * @param {Float} atomData.zCoord - z coordinate in Å
   * @param {?Float} atomData.bFactor - B-factor in in Å^2
   * @param {?Char} atomData.altLoc - alternate location identifier
   * @param {?Float} atomData.occupancy - occupancy of the atom
   */

  /**
   * @callback module:MmtfTraverse.onBond
   * @param {Object} bondData
   * @param {Integer} bondData.atomIndex1 - index of the first atom
   * @param {Integer} bondData.atomIndex2 - index of the secound atom
   * @param {Integer} bondData.bondOrder - bond order, allowed values are 1 to 3
   */


  /**
   * Traverse the MMTF structure data.
   * @static
   * @param {module:MmtfDecode.MmtfData} mmtfData - decoded mmtf data
   * @param {Object} eventCallbacks
   * @param {module:MmtfTraverse.onModel} [eventCallbacks.onModel] - called for each model
   * @param {module:MmtfTraverse.onChain} [eventCallbacks.onChain] - called for each chain
   * @param {module:MmtfTraverse.onGroup} [eventCallbacks.onGroup] - called for each group
   * @param {module:MmtfTraverse.onAtom} [eventCallbacks.onAtom] - called for each atom
   * @param {module:MmtfTraverse.onBond} [eventCallbacks.onBond] - called for each bond
   * @param {Object} [params] - traversal parameters
   * @param {Boolean} [params.firstModelOnly] - traverse only the first model
   */
  function traverseMmtf( mmtfData, eventCallbacks, params ){

      params = params || {};

      var firstModelOnly = params.firstModelOnly;

      // setup callbacks
      var onModel = eventCallbacks.onModel;
      var onChain = eventCallbacks.onChain;
      var onGroup = eventCallbacks.onGroup;
      var onAtom = eventCallbacks.onAtom;
      var onBond = eventCallbacks.onBond;

      // setup index counters
      var modelIndex = 0;
      var chainIndex = 0;
      var groupIndex = 0;
      var atomIndex = 0;

      var modelFirstAtomIndex = 0;
      var modelLastAtomIndex = -1;

      // setup optional fields
      var chainNameList = mmtfData.chainNameList;
      var secStructList = mmtfData.secStructList;
      var insCodeList = mmtfData.insCodeList;
      var sequenceIndexList = mmtfData.sequenceIndexList;
      var atomIdList = mmtfData.atomIdList;
      var bFactorList = mmtfData.bFactorList;
      var altLocList = mmtfData.altLocList;
      var occupancyList = mmtfData.occupancyList;
      var bondAtomList = mmtfData.bondAtomList;
      var bondOrderList = mmtfData.bondOrderList;

      // hoisted loop variables
      var o, ol, i, j, k, kl;

      // loop over all models
      for( o = 0, ol = mmtfData.chainsPerModel.length; o < ol; ++o ){

          if( firstModelOnly && modelIndex > 0 ) break;

          var modelChainCount = mmtfData.chainsPerModel[ modelIndex ];

          if( onModel ){
              onModel({
                  chainCount: modelChainCount,
                  modelIndex: modelIndex
              });
          }

          for( i = 0; i < modelChainCount; ++i ){

              var chainGroupCount = mmtfData.groupsPerChain[ chainIndex ];
              if( onChain ){
                  var chainId = fromCharCode(
                      mmtfData.chainIdList.subarray( chainIndex * 4, chainIndex * 4 + 4 )
                  );
                  var chainName = null;
                  if( chainNameList ){
                      chainName = fromCharCode(
                          chainNameList.subarray( chainIndex * 4, chainIndex * 4 + 4 )
                      );
                  }
                  onChain({
                      groupCount: chainGroupCount,
                      chainIndex: chainIndex,
                      modelIndex: modelIndex,
                      chainId: chainId,
                      chainName: chainName
                  });
              }

              for( j = 0; j < chainGroupCount; ++j ){

                  var groupData = mmtfData.groupList[ mmtfData.groupTypeList[ groupIndex ] ];
                  var groupAtomCount = groupData.atomNameList.length;
                  if( onGroup ){
                      var secStruct = null;
                      if( secStructList ){
                          secStruct = secStructList[ groupIndex ];
                      }
                      var insCode = null;
                      if( mmtfData.insCodeList ){
                          insCode = String.fromCharCode( insCodeList[ groupIndex ] );
                      }
                      var sequenceIndex = null;
                      if( sequenceIndexList ){
                          sequenceIndex = sequenceIndexList[ groupIndex ];
                      }
                      onGroup({
                          atomCount: groupAtomCount,
                          groupIndex: groupIndex,
                          chainIndex: chainIndex,
                          modelIndex: modelIndex,
                          groupId: mmtfData.groupIdList[ groupIndex ],
                          groupType: mmtfData.groupTypeList[ groupIndex ],
                          groupName: groupData.groupName,
                          singleLetterCode: groupData.singleLetterCode,
                          chemCompType: groupData.chemCompType,
                          secStruct: secStruct,
                          insCode: insCode,
                          sequenceIndex: sequenceIndex
                      });
                  }

                  for( k = 0; k < groupAtomCount; ++k ){

                      if( onAtom ){
                          var atomId = null;
                          if( atomIdList ){
                              atomId = atomIdList[ atomIndex ];
                          }
                          var bFactor = null;
                          if( bFactorList ){
                              bFactor = bFactorList[ atomIndex ];
                          }
                          var altLoc = null;
                          if( altLocList ){
                              altLoc = String.fromCharCode( altLocList[ atomIndex ] );
                          }
                          var occupancy = null;
                          if( occupancyList ){
                              occupancy = occupancyList[ atomIndex ];
                          }
                          onAtom({
                              atomIndex: atomIndex,
                              groupIndex: groupIndex,
                              chainIndex: chainIndex,
                              modelIndex: modelIndex,
                              atomId: atomId,
                              element: groupData.elementList[ k ],
                              atomName: groupData.atomNameList[ k ],
                              formalCharge: groupData.formalChargeList[ k ],
                              xCoord: mmtfData.xCoordList[ atomIndex ],
                              yCoord: mmtfData.yCoordList[ atomIndex ],
                              zCoord: mmtfData.zCoordList[ atomIndex ],
                              bFactor: bFactor,
                              altLoc: altLoc,
                              occupancy: occupancy
                          });
                      }

                      atomIndex += 1;
                  }

                  if( onBond ){
                      // intra group bonds
                      var groupBondAtomList = groupData.bondAtomList;
                      for( k = 0, kl = groupData.bondOrderList.length; k < kl; ++k ){
                          onBond({
                              atomIndex1: atomIndex - groupAtomCount + groupBondAtomList[ k * 2 ],
                              atomIndex2: atomIndex - groupAtomCount + groupBondAtomList[ k * 2 + 1 ],
                              bondOrder: groupData.bondOrderList[ k ]
                          });
                      }
                  }

                  groupIndex += 1;
              }

              chainIndex += 1;
          }

          modelFirstAtomIndex = modelLastAtomIndex + 1;
          modelLastAtomIndex = atomIndex - 1;  // subtract one as it already has been incremented

          if( onBond ){
              // inter group bonds
              if( bondAtomList ){
                  for( k = 0, kl = bondAtomList.length; k < kl; k += 2 ){
                      var atomIndex1 = bondAtomList[ k ];
                      var atomIndex2 = bondAtomList[ k + 1 ];
                      if( ( atomIndex1 >= modelFirstAtomIndex && atomIndex1 <= modelLastAtomIndex ) ||
                          ( atomIndex2 >= modelFirstAtomIndex && atomIndex2 <= modelLastAtomIndex )
                      ){
                          onBond({
                              atomIndex1: atomIndex1,
                              atomIndex2: atomIndex2,
                              bondOrder: bondOrderList ? bondOrderList[ k / 2 ] : null
                          });
                      }
                  }
              }
          }

          modelIndex += 1;
      }

  }

  /**
   * Version name
   * @static
   * @type {String}
   */
  var version = "v1.0.0";

  /**
   * Version name
   * @private
   * @type {String}
   */
  var baseUrl = "http://mmtf.rcsb.org/v1.0/";

  /**
   * URL of the RCSB webservice to obtain MMTF files
   * @static
   * @type {String}
   */
  var fetchUrl = baseUrl + "full/";

  /**
   * URL of the RCSB webservice to obtain reduced MMTF files
   * @static
   * @type {String}
   */
  var fetchReducedUrl = baseUrl + "reduced/";

  /**
   * Encode MMTF fields
   * @static
   * @param  {module:MmtfDecode.MmtfData} mmtfData - mmtf data
   * @return {Uint8Array} encoded MMTF fields
   */
  function encode( mmtfData ){
      return encodeMsgpack( encodeMmtf( mmtfData ) );
  }

  /**
   * Decode MMTF fields
   * @static
   * @example
   * // bin is Uint8Array containing the mmtf msgpack
   * var mmtfData = MMTF.decode( bin );
   * console.log( mmtfData.numAtoms );
   *
   * @param  {Uint8Array|ArrayBuffer|module:MmtfDecode.EncodedMmtfData} binOrDict - binary MessagePack or encoded MMTF data
   * @param  {Object} [params] - decoding parameters
   * @param {String[]} params.ignoreFields - names of optional fields not to decode
   * @return {module:MmtfDecode.MmtfData} mmtfData
   */
  function decode( binOrDict, params ){
  	// make sure binOrDict is not a plain Arraybuffer
      if( binOrDict instanceof ArrayBuffer ){
          binOrDict = new Uint8Array( binOrDict );
      }

      var inputDict;
      if( binOrDict instanceof Uint8Array ){
          // get dict from msgpack
          inputDict = decodeMsgpack( binOrDict );
      }else{
          // already a dict
          inputDict = binOrDict;
      }

      return decodeMmtf( inputDict, params );
  }

  /**
   * @callback module:MMTF.onLoad
   * @param {module:MmtfDecode.MmtfData} mmtfData - decoded mmtf data object
   */

  /**
   * helper method to fetch binary files from an URL
   * @private
   * @param  {String} pdbid - PDB ID to fetch
   * @param  {String} baseUrl - URL to fetch from
   * @param  {module:MMTF.onLoad} onLoad - callback( mmtfData )
   * @param  {Function} onError - callback( error )
   * @return {undefined}
   */
  function _fetch( pdbid, baseUrl, onLoad, onError ){
      var xhr = new XMLHttpRequest();
      function _onLoad(){
          try{
              var mmtfData = decode( xhr.response );
              onLoad( mmtfData );
          }catch( error ){
              onError( error );
          }
      }
      xhr.addEventListener( "load", _onLoad, true );
      xhr.addEventListener( "error", onError, true );
      xhr.responseType = "arraybuffer";
      xhr.open( "GET", baseUrl + pdbid.toUpperCase() );
      xhr.send();
  }

  /**
   * Fetch MMTF file from RCSB webservice which contains
   * @static
   * @example
   * MMTF.fetch(
   *     "3PQR",
   *     // onLoad callback
   *     function( mmtfData ){ console.log( mmtfData ) },
   *     // onError callback
   *     function( error ){ console.error( error ) }
   * );
   *
   * @param  {String} pdbid - PDB ID to fetch
   * @param  {module:MMTF.onLoad} onLoad - callback( mmtfData )
   * @param  {Function} onError - callback( error )
   * @return {undefined}
   */
  function fetch( pdbid, onLoad, onError ){
      _fetch( pdbid, fetchUrl, onLoad, onError );
  }

  /**
   * Fetch reduced MMTF file from RCSB webservice which contains
   * protein C-alpha, nucleotide phosphate and ligand atoms
   * @static
   * @example
   * MMTF.fetchReduced(
   *     "3PQR",
   *     // onLoad callback
   *     function( mmtfData ){ console.log( mmtfData ) },
   *     // onError callback
   *     function( error ){ console.error( error ) }
   * );
   *
   * @param  {String} pdbid - PDB ID to fetch
   * @param  {module:MMTF.onLoad} onLoad - callback( mmtfData )
   * @param  {Function} onError - callback( error )
   * @return {undefined}
   */
  function fetchReduced( pdbid, onLoad, onError ){
      _fetch( pdbid, fetchReducedUrl, onLoad, onError );
  }

  exports.encode = encode;
  exports.decode = decode;
  exports.traverse = traverseMmtf;
  exports.fetch = fetch;
  exports.fetchReduced = fetchReduced;
  exports.version = version;
  exports.fetchUrl = fetchUrl;
  exports.fetchReducedUrl = fetchReducedUrl;
  exports.encodeMsgpack = encodeMsgpack;
  exports.encodeMmtf = encodeMmtf;
  exports.decodeMsgpack = decodeMsgpack;
  exports.decodeMmtf = decodeMmtf;

  Object.defineProperty(exports, '__esModule', { value: true });

}));
